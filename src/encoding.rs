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
use crate::{RI, SPos, State, StateMap, Tables, X, XArr};

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
