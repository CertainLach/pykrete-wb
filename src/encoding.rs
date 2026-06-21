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

use std::io::{self, Read, Write};

use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use rand::seq::SliceRandom;

use crate::consts::INV_SHIFT_ROWS_TAB;
use crate::consts::SHIFT_ROWS_TAB;
use crate::internal::XorEncodingSingle;
use crate::mat::{GF2, MatGF2, VecGF2};
use crate::xor::Bijection4;
use crate::{
	FoldIdx, FoldMap, HighLow, HighLowMap, NibbleMap, RI, SPos, State, StateMap, Tables, U4, X,
	XArr,
};

#[derive(Clone, Copy)]
pub struct Bijection8(XArr<X>);
impl Bijection8 {
	fn identity() -> Self {
		Self(XArr::from_fn(|i| i))
	}
	/// map if inv == false, unmap if true
	fn mapunmap(&self, v: X, inv: bool) -> X {
		if inv { self.unmap(v) } else { self.map(v) }
	}
	fn map(&self, v: X) -> X {
		self.0[v]
	}
	fn unmap(&self, v: X) -> X {
		for x in X::all() {
			if self.map(x) == v {
				return x;
			}
		}
		unreachable!()
	}
}
impl Bijection8 {
	pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for ele in self.0.0.iter() {
			out.write_all(&[ele.0])?
		}
		Ok(())
	}
	pub fn read_from(mut input: impl Read) -> io::Result<Self> {
		Ok(Self(XArr::try_from_fn(|_| {
			let mut v = [0];
			input.read_exact(&mut v)?;
			Ok(X(v[0])) as io::Result<_>
		})?))
	}
}
impl Default for Bijection8 {
	fn default() -> Self {
		Self::identity()
	}
}

impl Distribution<Bijection8> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Bijection8 {
		let mut out = Bijection8::identity();
		out.0.0.shuffle(rng);
		out
	}
}

pub struct ExternalEncoding(StateMap<Bijection8>);
impl ExternalEncoding {
	pub fn identity() -> Self {
		Self(Default::default())
	}
	pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for ele in self.0.0.iter() {
			ele.write_to(&mut out)?;
		}
		Ok(())
	}
	pub fn read_from(mut input: impl Read) -> io::Result<Self> {
		Ok(Self(StateMap::try_from_fn(|_| {
			Bijection8::read_from(&mut input)
		})?))
	}
}

impl Distribution<ExternalEncoding> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ExternalEncoding {
		ExternalEncoding(rng.random())
	}
}

pub fn apply_encoding(state: &mut State, encoding: &ExternalEncoding, inv: bool) {
	for pos in SPos::all() {
		let original_value = state[pos];
		state[pos] = encoding.0[pos].mapunmap(original_value, inv);
	}
}

// TODO: Provide better aliases (decode-input, encode-output) back
impl<const NRM1: usize> Tables<NRM1> {
	pub fn input_encoding(&mut self, e: &ExternalEncoding, inv: bool) {
		let tybox = &mut self.tyboxes.0[RI(0)];
		let shift_tab = if self.inv {
			&INV_SHIFT_ROWS_TAB
		} else {
			&SHIFT_ROWS_TAB
		};

		for pos in SPos::all() {
			let table_copy = tybox.0[pos];
			let table = &mut tybox.0[pos];
			for encoded_x in X::all() {
				// cipher performs shift_rows as a first step, input decoding should account for that.
				let original_pos = shift_tab.map(pos);
				let x = e.0[original_pos].mapunmap(encoded_x, inv);
				let temp = table_copy[x];
				table[encoded_x] = temp;
			}
		}
	}

	pub fn output_encoding(&mut self, e: &ExternalEncoding, inv: bool) {
		let tbox = &mut self.tboxes_last;
		for pos in SPos::all() {
			let table_copy = XArr(X::ALL.map(|x| tbox.0[x][pos]));
			for x in X::all() {
				let temp = table_copy[x];
				let encoded = e.0[pos].mapunmap(temp, inv);
				tbox.0[x][pos] = encoded;
			}
		}
	}
}

pub struct LinearExternalEncoding(MatGF2<128>);
impl LinearExternalEncoding {
	pub fn identity() -> Self {
		Self(MatGF2::<128>::identity_n())
	}
	pub fn random<R: Rng>(rng: &mut R) -> Self {
		Self(MatGF2::<128>::random_invertible(rng))
	}
	fn matrix(&self, inv: bool) -> MatGF2<128> {
		if inv {
			self.0
				.inverse()
				.expect("encoding matrix must be invertible")
		} else {
			self.0
		}
	}
	pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for i in 0..128 {
			let mut row = [0u8; 16];
			for j in 0..128 {
				if self.0.get(i, j) == GF2::ONE {
					row[j / 8] |= 1 << (7 - (j % 8));
				}
			}
			out.write_all(&row)?;
		}
		Ok(())
	}
	pub fn read_from(mut input: impl Read) -> io::Result<Self> {
		let mut m = MatGF2::<128>::new();
		for i in 0..128 {
			let mut row = [0u8; 16];
			input.read_exact(&mut row)?;
			for j in 0..128 {
				if (row[j / 8] >> (7 - (j % 8))) & 1 == 1 {
					m.set(i, j, GF2::ONE);
				}
			}
		}
		Ok(Self(m))
	}
}
impl Default for LinearExternalEncoding {
	fn default() -> Self {
		Self::identity()
	}
}

fn state_to_vec(state: &State) -> VecGF2<128> {
	let mut v = VecGF2::<128>::new();
	for (idx, pos) in SPos::all().enumerate() {
		let b = state[pos].0;
		for bit in 0..8 {
			v.set(idx * 8 + bit, GF2((b >> (7 - bit)) & 1 == 1));
		}
	}
	v
}
fn vec_to_state(v: &VecGF2<128>) -> State {
	State::from_fn(|pos| {
		let mut b = 0u8;
		for bit in 0..8 {
			if v.get(pos.as_index() * 8 + bit) == GF2::ONE {
				b |= 1 << (7 - bit);
			}
		}
		X(b)
	})
}

fn x_nib(x: X, hl: HighLow) -> U4 {
	let (hi, lo) = x.as_nibs();
	match hl {
		HighLow::High => hi,
		HighLow::Low => lo,
	}
}

#[derive(Debug)]
pub struct LinearNetworkOutput(pub(crate) StateMap<XArr<X>>);
impl LinearNetworkOutput {
	pub(crate) fn apply(&self, state: &mut State) {
		*state = State::from_fn(|q| self.0[q][state[q]]);
	}
}

#[derive(Clone, Debug)]
pub struct LinearNetwork {
	pub(crate) in_tab: StateMap<XArr<State>>,
	pub(crate) fold: StateMap<HighLowMap<FoldMap<NibbleMap<Bijection4>>>>,
}
impl LinearNetwork {
	fn build<R: Rng>(
		mat: &MatGF2<128>,
		decode: impl Fn(SPos, X) -> X,
		root: &XorEncodingSingle,
		rng: &mut R,
	) -> Self {
		// Private operand encoding for the partial from input ip at output q.
		let pin = StateMap::from_fn(|_| {
			StateMap::from_fn(|_| HighLowMap::from_fn(|_| rng.random::<Bijection4>()))
		});

		let in_tab = StateMap::from_fn(|ip| {
			XArr::from_fn(|byte| {
				let c = decode(ip, byte);
				let mut single = State::default();
				single[ip] = c;
				let partial = vec_to_state(&mat.mul_vec(&state_to_vec(&single)));
				State::from_fn(|q| {
					let (rhi, rlo) = partial[q].as_nibs();
					let hi = pin[ip][q][HighLow::High].map(rhi);
					let lo = pin[ip][q][HighLow::Low].map(rlo);
					X::nibs(hi, lo)
				})
			})
		});

		let fold = StateMap::from_fn(|q| {
			HighLowMap::from_fn(|hl| {
				let root_enc = root.nibble(q, hl);
				let mut left = pin[SPos::from_index(0)][q][hl];
				FoldMap::from_fn(|k| {
					let ip = k.next_spos();
					let right = pin[ip][q][hl];
					let out_enc = if k.is_last() { root_enc } else { rng.random() };
					let out = NibbleMap::from_fn(|a| {
						Bijection4::new(
							U4::ALL.map(|b| out_enc.map(left.unmap(a) ^ right.unmap(b))),
						)
					});
					left = out_enc;
					out
				})
			})
		});

		Self { in_tab, fold }
	}

	pub(crate) fn apply(&self, state: &mut State) {
		let input = *state;
		*state = State::from_fn(|q| {
			let nib = HighLowMap::from_fn(|hl| {
				let ip0 = SPos::from_index(0);
				let mut acc = x_nib(self.in_tab[ip0][input[ip0]][q], hl);
				let tree = &self.fold[q][hl];
				for k in FoldIdx::all() {
					let ip = k.next_spos();
					let b = x_nib(self.in_tab[ip][input[ip]][q], hl);
					acc = tree[k][acc].map(b);
				}
				acc
			});
			X::nibs(*nib.high(), *nib.low())
		});
	}
}

impl<const NRM1: usize> Tables<NRM1> {
	pub fn input_linear_encoding<R: Rng>(
		&mut self,
		linear: &LinearExternalEncoding,
		nonlinear: Option<&ExternalEncoding>,
		inv: bool,
		rng: &mut R,
	) {
		assert!(
			self.input_linear.is_none(),
			"linear input encoding already installed"
		);
		let shift = if self.inv {
			INV_SHIFT_ROWS_TAB
		} else {
			SHIFT_ROWS_TAB
		};
		let mat = linear.matrix(inv);
		let boundary: XorEncodingSingle = rng.random();
		let net = LinearNetwork::build(
			&mat,
			|pos, b| match nonlinear {
				Some(n) => n.0[pos].unmap(b),
				None => b,
			},
			&boundary,
			rng,
		);
		self.tyboxes.0[RI(0)].encode(&Some(boundary), &None, &shift);
		self.input_linear = Some(net);
	}

	pub fn output_linear_encoding<R: Rng>(
		&mut self,
		linear: &LinearExternalEncoding,
		nonlinear: Option<&ExternalEncoding>,
		inv: bool,
		rng: &mut R,
	) {
		assert!(
			self.output_linear.is_none(),
			"linear output encoding already installed"
		);
		let mat = linear.matrix(inv);
		let last_enc: XorEncodingSingle = rng.random();
		for (_x, out) in self.tboxes_last.0.iter_mut() {
			for pos in SPos::all() {
				out[pos] = last_enc.map(pos, out[pos]);
			}
		}
		let (net, netout) = match nonlinear {
			Some(n) => {
				let boundary: XorEncodingSingle = rng.random();
				let net =
					LinearNetwork::build(&mat, |pos, b| last_enc.unmap(pos, b), &boundary, rng);
				(
					net,
					Some(LinearNetworkOutput(StateMap::from_fn(|pos| {
						XArr::from_fn(|b| n.0[pos].map(boundary.unmap(pos, b)))
					}))),
				)
			}
			None => (
				LinearNetwork::build(
					&mat,
					|pos, b| last_enc.unmap(pos, b),
					&XorEncodingSingle::identity(),
					rng,
				),
				None,
			),
		};
		self.output_linear = Some(net);
		self.output_linear_out = netout;
	}
}

pub fn apply_linear_encoding(state: &mut State, encoding: &LinearExternalEncoding, inv: bool) {
	let v = encoding.matrix(inv).mul_vec(&state_to_vec(state));
	*state = vec_to_state(&v);
}
