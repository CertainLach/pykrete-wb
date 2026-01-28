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

use tracing::{Level, debug, instrument};

use crate::consts::{INV_SHIFT_ROWS_TAB, SHIFT_ROWS_TAB};
use crate::key::RoundKeys;
use crate::sbox::{PrecomputedSBox, SBox};
use crate::{RI, RIArr, SPos, State, X, XArr, add_round_key, add_shifted_round_key, sub_bytes};

#[derive(Clone, Copy, Default, Debug)]
pub struct Tbox(pub XArr<State>);
impl Tbox {
	pub(crate) fn apply(&self, state: &mut State) {
		for pos in SPos::all() {
			state[pos] = self.0[state[pos]][pos];
		}
	}
}

pub struct Tboxes<const NRM1: usize> {
	tboxes: RIArr<Tbox, NRM1>,
	pub last: Tbox,
}
impl<const NRM1: usize> Tboxes<NRM1> {
	pub(crate) fn new(tboxes: RIArr<Tbox, NRM1>, last: Tbox) -> Self {
		Self { tboxes, last }
	}
	pub fn from_round_keys<S: SBox + Copy>(
		sbox: S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
	) -> Self {
		if inv {
			Self::from_round_keys_inv(&PrecomputedSBox::inversed(&sbox), round_keys)
		} else {
			Self::from_round_keys_forward(sbox, round_keys)
		}
	}

	#[instrument(level = Level::DEBUG, name = "tboxes_from_round_inv", skip(round_keys, inv_sbox))]
	fn from_round_keys_inv<S: SBox + Copy>(inv_sbox: S, round_keys: &RoundKeys<NRM1>) -> Self {
		debug!("last (as in original) round");
		let mut tboxes: RIArr<Tbox, NRM1> = Default::default();
		let mut last: XArr<State> = Default::default();

		let (a, b) = &round_keys.last;
		for (x, xbox) in last.iter_mut() {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, b, &INV_SHIFT_ROWS_TAB);
			sub_bytes(inv_sbox, &mut state);
			add_round_key(&mut state, a);

			*xbox = state;
		}
		tboxes.0[0] = Tbox(last);

		debug!("inner rounds");
		for (i, tbox) in tboxes.iter_mut().skip(1).rev() {
			for (x, xbox) in tbox.0.iter_mut() {
				let mut state = State::from_x(x);
				sub_bytes(inv_sbox, &mut state);
				add_round_key(&mut state, &round_keys.rounds[i]);

				*xbox = state;
			}
		}
		tboxes.0[1..].reverse();

		debug!("first and last round");
		let last = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			sub_bytes(inv_sbox, &mut state);
			add_round_key(&mut state, &round_keys.rounds[RI(0)]);

			state
		}));

		Self {
			tboxes,
			last: Tbox(last),
		}
	}
	#[instrument(level = Level::DEBUG, name = "tboxes_from_round", skip(round_keys, sbox))]
	fn from_round_keys_forward<S: SBox + Copy>(sbox: S, round_keys: &RoundKeys<NRM1>) -> Self {
		let mut tboxes: RIArr<Tbox, NRM1> = Default::default();

		debug!("base rounds");
		for (r, rtbox) in tboxes.iter_mut() {
			for (x, xbox) in rtbox.0.iter_mut() {
				let mut state = State::from_x(x);
				add_shifted_round_key(&mut state, &round_keys.rounds[r], &SHIFT_ROWS_TAB);
				sub_bytes(sbox, &mut state);

				*xbox = state;
			}
		}

		debug!("last round");

		let (a, b) = &round_keys.last;
		let last: XArr<State> = XArr(X::ALL.map(|x| {
			let mut state = State::from_x(x);
			add_shifted_round_key(&mut state, a, &SHIFT_ROWS_TAB);
			sub_bytes(sbox, &mut state);
			add_round_key(&mut state, b);

			state
		}));

		Self {
			tboxes,
			last: Tbox(last),
		}
	}
	pub fn get(&self, r: RI, x: X, ij: SPos) -> X {
		self.tboxes[r].0[x][ij]
	}
}
