// Pykrete Whitebox
// Copyright (C) 2026  Yaroslav Bolyukin <iam@0la.ch>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License v3 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use rand::Rng;

use crate::dual::{Dual, Q};
use crate::key::{Key, RoundKeys};
use crate::mat::MatGF2;
use crate::sbox::SBox;
use crate::tbox::{Tbox, Tboxes};
use crate::ty::Ty;
use crate::tybox::{Work, WorkRounds};
use crate::xor::ShiftRowsBijection;
use crate::{RI, RIArr, RowMap, SPos, SRow, State, StateMap, Tables, X, XArr, add_round_key};

use crate::consts::{INV_SHIFT_ROWS_TAB, SHIFT_ROWS_TAB};
use crate::sbox::PrecomputedSBox;

#[derive(Debug)]
pub struct KarroumiConfig<const NRM1: usize> {
	pub rounds: RIArr<Dual, NRM1>,
	pub last: Dual,
}

impl<const NRM1: usize> KarroumiConfig<NRM1> {
	pub fn new(rounds: RIArr<Dual, NRM1>, last: Dual) -> Self {
		Self { rounds, last }
	}

	pub fn random<R: Rng>(rng: &mut R) -> Self {
		let rounds = RIArr::from_fn(|_| Dual::random(rng));

		let last = { Dual::random(rng) };

		Self { rounds, last }
	}

	pub fn standard() -> Self {
		Self {
			rounds: RIArr::from_fn(|_| Dual::STANDARD),
			last: Dual::STANDARD,
		}
	}

	pub fn round_config(&self, r: RI) -> Dual {
		if r.0 < NRM1 {
			self.rounds[r]
		} else {
			self.last
		}
	}
}

fn expand_karroumi_keys<S: SBox, const NK: usize, const NRM1: usize>(
	sbox_fn: &impl Fn(Dual) -> S,
	key: Key<NK>,
	config: &KarroumiConfig<NRM1>,
) -> RoundKeys<NRM1> {
	use crate::dual::compute_rcon;
	use crate::{Word, rot_word, sub_word};

	const NB: usize = 4;
	let mut out = RoundKeys::new();

	let r0_q = Q::for_dual(config.rounds[RI(0)]);
	let mut transformed_key = key;
	transformed_key.apply_q(r0_q);

	for i in 0..NK {
		out.set_schedule_word(i, transformed_key.word(i));
	}

	let mut rcon_idx = 1usize;
	for i in NK..(NB * (NRM1 + 2)) {
		let w_im1 = out.get_schedule_word(i - 1);
		let w_imn = out.get_schedule_word(i - NK);

		let round_idx = i / NB;
		let round_config = if round_idx == 0 {
			config.rounds[RI(0)]
		} else if round_idx <= NRM1 {
			config.round_config(RI(round_idx - 1))
		} else {
			config.last
		};

		let sbox = sbox_fn(round_config);

		let w_i = w_imn
			^ if i.is_multiple_of(NK) {
				let rc = compute_rcon(rcon_idx, round_config);
				let r = sub_word(&sbox, rot_word(w_im1)) ^ Word::from_bytes([rc, 0, 0, 0]);
				rcon_idx += 1;
				r
			} else if NK > 6 && (i % NK) == 4 {
				sub_word(&sbox, w_im1)
			} else {
				w_im1
			};

		out.set_schedule_word(i, w_i);
	}

	out
}

pub struct KarroumiTy<const NRM1: usize>(RIArr<Ty, NRM1>);

impl<const NRM1: usize> KarroumiTy<NRM1> {
	pub fn new(inv: bool, config: &KarroumiConfig<NRM1>) -> Self {
		Self(if inv {
			RIArr::from_fn(|dec_idx| {
				let dual = if dec_idx.0 == 0 {
					config.last
				} else {
					config.rounds[RI(NRM1 - dec_idx.0)]
				};
				Ty::new(inv, dual)
			})
		} else {
			RIArr::from_fn(|r| Ty::new(inv, config.rounds[r]))
		})
	}

	pub fn get(&self, r: RI) -> &Ty {
		&self.0[r]
	}
}

#[derive(Clone, Copy)]
pub struct Delta(MatGF2<8>);
impl Delta {
	pub fn new(mat: MatGF2<8>) -> Self {
		Self(mat)
	}
	pub fn apply(&self, x: X) -> X {
		X(self.0.apply(x.0))
	}
}

impl Tbox {
	fn apply_delta(&mut self, delta: Delta) {
		let old = self.0;
		for x in X::all() {
			self.0[x] = old[delta.apply(x)];
		}
	}
	fn apply_delta_inv(&mut self, delta: Delta) {
		for x in X::all() {
			for pos in SPos::all() {
				let old = self.0[x][pos];
				self.0[x][pos] = delta.apply(old);
			}
		}
	}
}
impl Work {
	fn apply_delta(&mut self, delta: Delta) {
		for pos in SPos::all() {
			let old = self.0[pos];
			self.0[pos] = XArr::from_fn(|x| old[delta.apply(x)]);
		}
	}
}

impl<const NRM1: usize> Tboxes<NRM1> {
	pub fn from_karroumi_round_keys<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
		config: &KarroumiConfig<NRM1>,
	) -> Self {
		if inv {
			Self::from_karroumi_round_keys_inv(sbox_fn, round_keys, config)
		} else {
			Self::from_karroumi_round_keys_forward(sbox_fn, round_keys, config)
		}
	}

	fn from_karroumi_round_keys_forward<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		config: &KarroumiConfig<NRM1>,
	) -> Self {
		use crate::{add_shifted_round_key, sub_bytes};

		let tboxes = RIArr(std::array::from_fn(|r| {
			let r = RI(r);
			let sbox = sbox_fn(config.rounds[r]);
			let mut rtbox = Tbox::default();
			for (x, xbox) in rtbox.0.iter_mut() {
				let mut state = State::from_x(x);
				add_shifted_round_key(&mut state, &round_keys.rounds[r], &SHIFT_ROWS_TAB);
				sub_bytes(&sbox, &mut state);
				*xbox = state;
			}
			rtbox
		}));

		let (a, b) = &round_keys.last;
		let sbox = sbox_fn(config.last);
		let last: XArr<State> = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, a, &SHIFT_ROWS_TAB);
			sub_bytes(&sbox, &mut state);
			add_round_key(&mut state, b);
			state
		}));

		Self::new(tboxes, Tbox(last))
	}

	fn from_karroumi_round_keys_inv<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		config: &KarroumiConfig<NRM1>,
	) -> Self {
		use crate::{add_round_key, add_shifted_round_key, sub_bytes};

		let mut tboxes_arr: [Tbox; NRM1] = std::array::from_fn(|_| Tbox::default());

		let (a, b) = &round_keys.last;
		let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(config.last));
		let mut last_as_first: XArr<State> = Default::default();
		for (x, xbox) in last_as_first.iter_mut() {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, b, &INV_SHIFT_ROWS_TAB);
			sub_bytes(&inv_sbox, &mut state);
			add_round_key(&mut state, a);
			*xbox = state;
		}
		tboxes_arr[0] = Tbox(last_as_first);

		for dec_idx in 1..NRM1 {
			let enc_round = NRM1 - dec_idx; // Encryption round being inverted
			let round_config = config.rounds[RI(enc_round)];
			let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(round_config));
			for (x, xbox) in tboxes_arr[dec_idx].0.iter_mut() {
				let mut state = State::from_x(x);
				sub_bytes(&inv_sbox, &mut state);
				add_round_key(&mut state, &round_keys.rounds[RI(enc_round)]);
				*xbox = state;
			}
		}

		let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(config.rounds[RI(0)]));
		let last = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			sub_bytes(&inv_sbox, &mut state);
			add_round_key(&mut state, &round_keys.rounds[RI(0)]);
			state
		}));

		Self::new(RIArr(tboxes_arr), Tbox(last))
	}
}

impl<const NRM1: usize> WorkRounds<NRM1> {
	pub fn new_karroumi_tyi(tboxes: &Tboxes<NRM1>, ty: &KarroumiTy<NRM1>) -> Self {
		Self(RIArr(std::array::from_fn(|r| {
			let r = RI(r);
			Work(StateMap::from_fn(|pos| {
				XArr::from_fn(|x| {
					let i = pos.column();
					let tboxv = tboxes.get(r, x, pos);
					ty.get(r).get_column(tboxv, i)
				})
			}))
		})))
	}
}

impl<const NRM1: usize> Tables<NRM1> {
	pub fn from_karroumi_key<S: SBox, const NK: usize>(
		sbox_fn: impl Fn(Dual) -> S,
		key: Key<NK>,
		inv: bool,
		base: Dual,
		config: &KarroumiConfig<NRM1>,
	) -> Self {
		let round_keys = expand_karroumi_keys(&sbox_fn, key, config);
		Self::from_karroumi_round_keys(sbox_fn, &round_keys, inv, base, config)
	}

	pub fn from_karroumi_round_keys<S: SBox>(
		sbox_fn: impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
		base: Dual,
		config: &KarroumiConfig<NRM1>,
	) -> Self {
		let tboxes = Tboxes::from_karroumi_round_keys(&sbox_fn, round_keys, inv, config);
		let ty = KarroumiTy::new(inv, config);

		let mut tyboxes = WorkRounds::new_karroumi_tyi(&tboxes, &ty);

		let mut tboxes_last = tboxes.last;

		let standard_q = Q::for_dual(Dual::STANDARD);
		let base_q = Q::for_dual(base);
		let last_q = Q::for_dual(config.last);
		let r0_q = Q::for_dual(config.rounds[RI(0)]);
		if inv {
			tyboxes.0[RI(0)].apply_delta(base_q.delta(last_q));

			for dec_idx in 1..NRM1 {
				let prev_repr = if dec_idx == 1 {
					config.last
				} else {
					config.rounds[RI(NRM1 - dec_idx + 1)]
				};
				tyboxes.0[RI(dec_idx)]
					.apply_delta(prev_repr.delta(config.rounds[RI(NRM1 - dec_idx)]));
			}

			tboxes_last.apply_delta(Q::for_dual(config.rounds[RI(1)]).delta(r0_q));

			tboxes_last.apply_delta_inv(r0_q.delta(standard_q));
		} else {
			{
				tyboxes.0[RI(0)].apply_delta(standard_q.delta(r0_q));
			}

			for r in 1..NRM1 {
				tyboxes.0[RI(r)].apply_delta(config.rounds[RI(r - 1)].delta(config.rounds[RI(r)]));
			}

			tboxes_last.apply_delta(Q::for_dual(config.rounds[RI(NRM1 - 1)]).delta(last_q));

			tboxes_last.apply_delta_inv(last_q.delta(base_q));
		}

		Self::new_base(tyboxes, tboxes_last, inv)
	}
}

// Can be replaced with RowMap<Delta>, but then it would be necessary
// to perform computations with ShiftRows in mind, and it would be more complicated
struct Delta4(StateMap<Delta>);

#[derive(Debug, Clone, Copy)]
pub struct Dual4(RowMap<Dual>);
impl Dual4 {
	const STANDARD: Self = Self::uniform(Dual::STANDARD);
	fn random<R: Rng>(rng: &mut R) -> Self {
		Self(RowMap::from_fn(|_| Dual::random(rng)))
	}
	pub const fn uniform(dual: Dual) -> Self {
		Self(RowMap([dual, dual, dual, dual]))
	}

	fn delta(&self, curr: Dual4, shift_rows: &ShiftRowsBijection) -> Delta4 {
		Delta4(StateMap::from_fn(|pos| {
			let row = pos.row();
			let prev_row = shift_rows.map(pos).row();

			self.0[prev_row].delta(curr.0[row])
		}))
	}

	fn delta_inv(&self, to: Dual4) -> Delta4 {
		Delta4(StateMap::from_fn(|pos| {
			let row = pos.row();
			self.0[row].delta(to.0[row])
		}))
	}
}

#[derive(Debug)]
pub struct KarroumiConfig4<const NRM1: usize> {
	pub rounds: RIArr<Dual4, NRM1>,
	pub last: Dual4,
}

impl<const NRM1: usize> KarroumiConfig4<NRM1> {
	pub fn new(rounds: RIArr<Dual4, NRM1>, last: Dual4) -> Self {
		Self { rounds, last }
	}

	pub fn random<R: Rng>(rng: &mut R) -> Self {
		let rounds = RIArr::from_fn(|_| Dual4::random(rng));
		let last = Dual4::random(rng);
		Self { rounds, last }
	}

	pub fn standard() -> Self {
		Self {
			rounds: RIArr::from_fn(|_| Dual4::STANDARD),
			last: Dual4::STANDARD,
		}
	}

	pub fn round_config(&self, r: RI) -> Dual4 {
		if r.0 < NRM1 {
			self.rounds[r]
		} else {
			self.last
		}
	}
}

impl Tbox {
	fn apply_delta4(&mut self, deltas: Delta4) {
		let old = self.0;
		for x in X::all() {
			for pos in SPos::all() {
				let delta = deltas.0[pos];
				self.0[x][pos] = old[delta.apply(x)][pos];
			}
		}
	}

	fn apply_delta4_inv(&mut self, deltas: Delta4) {
		for x in X::all() {
			for pos in SPos::all() {
				let delta = deltas.0[pos];
				let old = self.0[x][pos];
				self.0[x][pos] = delta.apply(old);
			}
		}
	}
}

impl Work {
	fn apply_delta4(&mut self, deltas: Delta4) {
		for pos in SPos::all() {
			let delta = deltas.0[pos];
			let old = self.0[pos];
			self.0[pos] = XArr::from_fn(|x| old[X(delta.0.apply(x.0))]);
		}
	}
}

pub struct KarroumiTy4<const NRM1: usize>(RIArr<RowMap<Ty>, NRM1>);

impl<const NRM1: usize> KarroumiTy4<NRM1> {
	pub fn new(inv: bool, config: &KarroumiConfig4<NRM1>) -> Self {
		Self(if inv {
			RIArr::from_fn(|dec_idx| {
				let enc_config = if dec_idx.0 == 0 {
					&config.last
				} else {
					&config.rounds[RI(NRM1 - dec_idx.0)]
				};
				RowMap::from_fn(|row| Ty::new(inv, enc_config.0[row]))
			})
		} else {
			RIArr::from_fn(|r| RowMap::from_fn(|row| Ty::new(inv, config.rounds[r].0[row])))
		})
	}

	pub fn get(&self, r: RI, row: SRow) -> &Ty {
		&self.0[r][row]
	}
}

fn expand_karroumi_keys_4<const NK: usize, const NRM1: usize>(key: Key<NK>) -> RoundKeys<NRM1> {
	use crate::key::expand_nonstandard_keys;
	expand_nonstandard_keys(&PrecomputedSBox::aes_standard(), key, Dual::STANDARD)
}

impl<const NRM1: usize> Tboxes<NRM1> {
	pub fn from_karroumi4_round_keys<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
		config: &KarroumiConfig4<NRM1>,
	) -> Self {
		if inv {
			Self::from_karroumi4_round_keys_inv(sbox_fn, round_keys, config)
		} else {
			Self::from_karroumi4_round_keys_forward(sbox_fn, round_keys, config)
		}
	}

	fn from_karroumi4_round_keys_forward<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		config: &KarroumiConfig4<NRM1>,
	) -> Self {
		use crate::{add_round_key, add_shifted_round_key};

		let tboxes = RIArr(std::array::from_fn(|r| {
			let r = RI(r);
			let mut rtbox = Tbox::default();

			let mut transformed_rk = round_keys.rounds[r].clone();
			for key_pos in SPos::all() {
				let state_pos = INV_SHIFT_ROWS_TAB.map(key_pos);
				let q = Q::for_dual(config.rounds[r].0[state_pos.row()]);
				transformed_rk.0[key_pos.0] = q.apply(round_keys.rounds[r].0[key_pos.0]);
			}

			for (x, xbox) in rtbox.0.iter_mut() {
				let mut state = State::from_x(x);
				add_shifted_round_key(&mut state, &transformed_rk, &SHIFT_ROWS_TAB);
				for pos in SPos::all() {
					let sbox = sbox_fn(config.rounds[r].0[pos.row()]);
					state[pos] = sbox.sub_byte(state[pos]);
				}
				*xbox = state;
			}
			rtbox
		}));

		let (a, b) = &round_keys.last;
		let mut transformed_a = a.clone();
		for key_pos in SPos::all() {
			let state_pos = INV_SHIFT_ROWS_TAB.map(key_pos);
			let q = Q::for_dual(config.last.0[state_pos.row()]);
			transformed_a.0[key_pos.0] = q.apply(a.0[key_pos.0]);
		}
		let mut transformed_b = b.clone();
		for pos in SPos::all() {
			let q = Q::for_dual(config.last.0[pos.row()]);
			transformed_b.0[pos.0] = q.apply(b.0[pos.0]);
		}

		let last: XArr<State> = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, &transformed_a, &SHIFT_ROWS_TAB);
			for pos in SPos::all() {
				let sbox = sbox_fn(config.last.0[pos.row()]);
				state[pos] = sbox.sub_byte(state[pos]);
			}
			add_round_key(&mut state, &transformed_b);
			state
		}));

		Self::new(tboxes, Tbox(last))
	}

	fn from_karroumi4_round_keys_inv<S: SBox>(
		sbox_fn: &impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		config: &KarroumiConfig4<NRM1>,
	) -> Self {
		use crate::{add_round_key, add_shifted_round_key};

		let mut tboxes_arr: [Tbox; NRM1] = std::array::from_fn(|_| Tbox::default());

		let (a, b) = &round_keys.last;

		let mut transformed_b = b.clone();
		for key_pos in SPos::all() {
			let state_pos = SHIFT_ROWS_TAB.map(key_pos);
			let q = Q::for_dual(config.last.0[state_pos.row()]);
			transformed_b.0[key_pos.0] = q.apply(b.0[key_pos.0]);
		}
		let mut transformed_a = a.clone();
		for pos in SPos::all() {
			let q = Q::for_dual(config.last.0[pos.row()]);
			transformed_a.0[pos.0] = q.apply(a.0[pos.0]);
		}

		let mut last_as_first: XArr<State> = Default::default();
		for (x, xbox) in last_as_first.iter_mut() {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, &transformed_b, &INV_SHIFT_ROWS_TAB);
			for pos in SPos::all() {
				let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(config.last.0[pos.row()]));
				state[pos] = inv_sbox.sub_byte(state[pos]);
			}
			add_round_key(&mut state, &transformed_a);
			*xbox = state;
		}
		tboxes_arr[0] = Tbox(last_as_first);

		for dec_idx in 1..NRM1 {
			let enc_round = NRM1 - dec_idx;
			let round_cfg = &config.rounds[RI(enc_round)];

			let mut transformed_rk = round_keys.rounds[RI(enc_round)].clone();
			for pos in SPos::all() {
				let q = Q::for_dual(round_cfg.0[pos.row()]);
				transformed_rk.0[pos.0] = q.apply(round_keys.rounds[RI(enc_round)].0[pos.0]);
			}

			for (x, xbox) in tboxes_arr[dec_idx].0.iter_mut() {
				let mut state = State::from_x(x);
				for pos in SPos::all() {
					let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(round_cfg.0[pos.row()]));
					state[pos] = inv_sbox.sub_byte(state[pos]);
				}
				add_round_key(&mut state, &transformed_rk);
				*xbox = state;
			}
		}

		let round_cfg = &config.rounds[RI(0)];
		let mut transformed_rk = round_keys.rounds[RI(0)].clone();
		for pos in SPos::all() {
			let q = Q::for_dual(round_cfg.0[pos.row()]);
			transformed_rk.0[pos.0] = q.apply(round_keys.rounds[RI(0)].0[pos.0]);
		}

		let last = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			for pos in SPos::all() {
				let inv_sbox = PrecomputedSBox::inversed(&sbox_fn(round_cfg.0[pos.row()]));
				state[pos] = inv_sbox.sub_byte(state[pos]);
			}
			add_round_key(&mut state, &transformed_rk);
			state
		}));

		Self::new(RIArr(tboxes_arr), Tbox(last))
	}
}

impl<const NRM1: usize> WorkRounds<NRM1> {
	pub fn new_karruomi4_tyi(tboxes: &Tboxes<NRM1>, ty: &KarroumiTy4<NRM1>) -> Self {
		Self(RIArr(std::array::from_fn(|r| {
			let r = RI(r);
			Work(StateMap::from_fn(|pos| {
				XArr::from_fn(|x| {
					let col = pos.column();
					let tboxv = tboxes.get(r, x, pos);
					ty.get(r, pos.row()).get_column(tboxv, col)
				})
			}))
		})))
	}
}

impl<const NRM1: usize> Tables<NRM1> {
	pub fn from_karroumi4_key<S: SBox, const NK: usize>(
		sbox_fn: impl Fn(Dual) -> S,
		key: Key<NK>,
		inv: bool,
		base: Dual4,
		config: &KarroumiConfig4<NRM1>,
	) -> Self {
		let round_keys = expand_karroumi_keys_4(key);
		Self::from_karroumi4_round_keys(sbox_fn, &round_keys, inv, base, config)
	}

	pub fn from_karroumi4_round_keys<S: SBox>(
		sbox_fn: impl Fn(Dual) -> S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
		base: Dual4,
		config: &KarroumiConfig4<NRM1>,
	) -> Self {
		let tboxes = Tboxes::from_karroumi4_round_keys(&sbox_fn, round_keys, inv, config);
		let ty = KarroumiTy4::new(inv, config);

		let mut tyboxes = WorkRounds::new_karruomi4_tyi(&tboxes, &ty);

		let mut tboxes_last = tboxes.last;

		let standard_configs = Dual4::STANDARD;

		if inv {
			let shift_rows = &INV_SHIFT_ROWS_TAB;
			tyboxes.0[RI(0)].apply_delta4(base.delta(config.last, shift_rows));

			for dec_idx in 1..NRM1 {
				let prev_repr = if dec_idx == 1 {
					config.last
				} else {
					config.rounds[RI(NRM1 - dec_idx + 1)]
				};
				tyboxes.0[RI(dec_idx)]
					.apply_delta4(prev_repr.delta(config.rounds[RI(NRM1 - dec_idx)], shift_rows));
			}

			tboxes_last.apply_delta4(config.rounds[RI(1)].delta(config.rounds[RI(0)], shift_rows));

			tboxes_last.apply_delta4_inv(config.rounds[RI(0)].delta_inv(standard_configs));
		} else {
			let shift_rows = &SHIFT_ROWS_TAB;
			tyboxes.0[RI(0)].apply_delta4(standard_configs.delta(config.rounds[RI(0)], shift_rows));

			for r in 1..NRM1 {
				tyboxes.0[RI(r)].apply_delta4(
					config.rounds[RI(r - 1)].delta(config.rounds[RI(r)], &SHIFT_ROWS_TAB),
				);
			}

			tboxes_last
				.apply_delta4(config.rounds[RI(NRM1 - 1)].delta(config.last, &SHIFT_ROWS_TAB));

			tboxes_last.apply_delta4_inv(config.last.delta_inv(base));
		}

		Self::new_base(WorkRounds(tyboxes.0), tboxes_last, inv)
	}
}

#[cfg(test)]
mod tests {
	use rand::rng;
	use test_case::test_case;

	use super::*;
	use crate::dual::IRREDUCIBLE_POLYNOMIALS;
	use crate::key::Aes128Key;

	fn make_sbox(config: Dual) -> PrecomputedSBox {
		PrecomputedSBox::for_dual(config)
	}

	#[test]
	fn standard_matches_original() {
		let base = Dual::STANDARD;
		let config = KarroumiConfig::<9>::standard();

		let tables = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
	}

	#[test]
	fn roundtrip() {
		let rng = &mut rng();

		let base = Dual::STANDARD;
		let config = KarroumiConfig::<9>::random(rng);

		let tables_enc = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let tables_dec = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		assert_ne!(data, original, "encrypted data should differ from original");

		tables_dec.cipher(&mut data);
		assert_eq!(data, original, "roundtrip should return original");
	}

	#[test]
	fn with_nonstandard_base() {
		let rng = &mut rng();

		let base = Dual::new(IRREDUCIBLE_POLYNOMIALS[1], 2);
		let config = KarroumiConfig::<9>::random(rng);

		let tables_enc = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let tables_dec = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);

		// Standard input (no Q applied before) - encryption produces base output
		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		assert_ne!(data, original, "encrypted should differ from original");

		// Decryption takes base input and produces standard output
		tables_dec.cipher(&mut data);
		assert_eq!(data, original, "roundtrip should return original");
	}

	#[test]
	fn nonstandard_base_standard_config_matches_aes() {
		let base = Dual::new(IRREDUCIBLE_POLYNOMIALS[1], 2);
		let config = KarroumiConfig::<9>::standard();

		let tables_enc = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;

		tables_enc.cipher(&mut data);

		let q_base = Q::for_dual(base);
		data.apply_q_inv(q_base);

		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
	}

	#[test]
	fn equivalence_to_standard_aes() {
		let rng = &mut rng();

		let base = Dual::STANDARD;
		let config = KarroumiConfig::<9>::random(rng);

		let karroumi_enc = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let karroumi_dec = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		karroumi_enc.cipher(&mut data);
		karroumi_dec.cipher(&mut data);

		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR);
	}

	#[test]
	fn full_security() {
		let rng = &mut rng();

		let base = Dual::STANDARD;
		let config = KarroumiConfig::<9>::random(rng);

		let tables_enc = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let tables_dec = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		tables_dec.cipher(&mut data);

		assert_eq!(data, original);
	}

	#[test]
	fn karroumi_aes256() {
		use crate::key::Aes256Key;

		let rng = &mut rng();

		let mut key = [0u8; 32];
		key[0] = 0x80;

		let base = Dual::STANDARD;
		let config = KarroumiConfig::<13>::random(rng);

		let tables_enc: Tables<13> =
			Tables::from_karroumi_key(make_sbox, Aes256Key::new(key), false, base, &config);
		let tables_dec: Tables<13> =
			Tables::from_karroumi_key(make_sbox, Aes256Key::new(key), true, base, &config);

		let original = State::from_bytes([0; 16]);
		let mut data = original;

		tables_enc.cipher(&mut data);
		tables_dec.cipher(&mut data);

		assert_eq!(data, original, "AES-256 Karroumi roundtrip should work");
	}

	// 4 Duals Per Round Tests

	#[test_case(false; "no security")]
	#[test_case(true; "security")]
	fn d4_standard_matches_original(security: bool) {
		let base = Dual4::STANDARD;
		let config = KarroumiConfig4::<9>::standard();

		let mut tables = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);

		if security {
			tables.apply_security(Security::full(), &mut rng())
		}

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
	}

	#[test_case(false; "no security")]
	#[test_case(true; "security")]
	fn d4_roundtrip(security: bool) {
		let rng = &mut rng();

		let base = Dual4::STANDARD;
		let config = KarroumiConfig4::<9>::random(rng);

		let mut tables_enc = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let mut tables_dec = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);
		if security {
			tables_enc.apply_security(Security::full(), rng);
			tables_dec.apply_security(Security::full(), rng);
		}

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		assert_ne!(data, original, "encrypted data should differ from original");

		tables_dec.cipher(&mut data);
		assert_eq!(data, original, "roundtrip should return original");
	}

	#[test_case(false; "no security")]
	#[test_case(true; "security")]
	fn d4_with_nonstandard_base(security: bool) {
		let rng = &mut rng();

		let base = Dual4(RowMap::from_fn(|col| {
			Dual::new(
				IRREDUCIBLE_POLYNOMIALS[col.as_index() + 1],
				col.as_index() % 8,
			)
		}));
		let config = KarroumiConfig4::<9>::random(rng);

		let mut tables_enc = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let mut tables_dec = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);
		if security {
			tables_enc.apply_security(Security::full(), rng);
			tables_dec.apply_security(Security::full(), rng);
		}

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		assert_ne!(data, original, "encrypted should differ from original");

		tables_dec.cipher(&mut data);
		assert_eq!(data, original, "roundtrip should return original");
	}

	#[test_case(false; "no security")]
	#[test_case(true; "security")]
	fn d4_aes256(security: bool) {
		use crate::key::Aes256Key;

		let rng = &mut rng();

		let mut key = [0u8; 32];
		key[0] = 0x80;

		let base = Dual4::STANDARD;
		let config = KarroumiConfig4::<13>::random(rng);

		let mut tables_enc: Tables<13> =
			Tables::from_karroumi4_key(make_sbox, Aes256Key::new(key), false, base, &config);
		let mut tables_dec: Tables<13> =
			Tables::from_karroumi4_key(make_sbox, Aes256Key::new(key), true, base, &config);
		if security {
			tables_enc.apply_security(Security::full(), rng);
			tables_dec.apply_security(Security::full(), rng);
		}

		let original = State::from_bytes([0; 16]);
		let mut data = original;

		tables_enc.cipher(&mut data);
		tables_dec.cipher(&mut data);

		assert_eq!(data, original, "AES-256 Karroumi4 roundtrip should work");
	}

	#[test_case(false; "no security")]
	#[test_case(true; "security")]
	fn d4_uniform_config_matches_karroumi1(security: bool) {
		let rng = &mut rng();

		let single_config = Dual::random(rng);
		let base = Dual4::STANDARD;

		let config1 = KarroumiConfig::<9> {
			rounds: RIArr::from_fn(|_| single_config),
			last: single_config,
		};
		let config4 = KarroumiConfig4::<9> {
			rounds: RIArr::from_fn(|_| Dual4::uniform(single_config)),
			last: Dual4::uniform(single_config),
		};

		let tables1 = Tables::from_karroumi_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			Dual::STANDARD,
			&config1,
		);
		let tables4 = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config4,
		);

		let mut data1 = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data4 = State::TWO_ONE_NINE_TWO_TEST_VECTOR;

		tables1.cipher(&mut data1);
		tables4.cipher(&mut data4);

		assert_eq!(data1, data4,);
	}

	#[test]
	fn d4_different_columns_per_round() {
		let config = KarroumiConfig4::<9> {
			rounds: RIArr::from_fn(|r| {
				Dual4(RowMap::from_fn(|col| {
					let poly_idx = (r.0 + col.as_index()) % 30;
					let power = (r.0 * 2 + col.as_index()) % 8;
					Dual::new(IRREDUCIBLE_POLYNOMIALS[poly_idx], power)
				}))
			}),
			last: Dual4(RowMap::from_fn(|col| {
				Dual::new(
					IRREDUCIBLE_POLYNOMIALS[col.as_index() + 10],
					col.as_index() % 8,
				)
			})),
		};

		let base = Dual4::STANDARD;

		let tables_enc = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			false,
			base,
			&config,
		);
		let tables_dec = Tables::from_karroumi4_key(
			make_sbox,
			Aes128Key::KUNG_FU_TEST_VECTOR,
			true,
			base,
			&config,
		);

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		tables_enc.cipher(&mut data);
		assert_ne!(data, original);

		tables_dec.cipher(&mut data);
		assert_eq!(data, original);
	}
}
