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
use rand::distr::{Distribution, StandardUniform};

use crate::tbox::Tbox;
use crate::tybox::Work;
use crate::xor::{Bijection4, ShiftRowsBijection, XorRound};
use crate::{
	ColumnMap, HighLow, HighLowMap, NibbleMap, Purpose, PurposeMap, RowMap, SColumn, SPos, SRow,
	U4, X,
};

#[derive(Clone)]
pub struct XorEncodingSingle(ColumnMap<ColumnMap<HighLowMap<Bijection4>>>);
impl XorEncodingSingle {
	pub fn identity() -> Self {
		Self(Default::default())
	}
	pub(crate) fn nibble(&self, pos: SPos, hl: HighLow) -> Bijection4 {
		self.0[pos.column()][SColumn(pos.row().0)][hl]
	}
	pub(crate) fn mapunmap(&self, pos: SPos, x: X, inv: bool) -> X {
		let (a, b) = x.as_nibs();
		let a = self.nibble(pos, HighLow::High).mapunmap(a, inv);
		let b = self.nibble(pos, HighLow::Low).mapunmap(b, inv);
		X::nibs(a, b)
	}
	pub(crate) fn map(&self, pos: SPos, x: X) -> X {
		self.mapunmap(pos, x, false)
	}
	pub(crate) fn unmap(&self, pos: SPos, x: X) -> X {
		self.mapunmap(pos, x, true)
	}
}
impl Default for XorEncodingSingle {
	fn default() -> Self {
		Self::identity()
	}
}
impl Distribution<XorEncodingSingle> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> XorEncodingSingle {
		let out: ColumnMap<ColumnMap<HighLowMap<Bijection4>>> = rng.random();
		XorEncodingSingle(out)
	}
}

pub struct XorEncoding {
	xor_high: PurposeMap<XorEncodingSingle>,
	xor_low: PurposeMap<XorEncodingSingle>,
	xor_output: PurposeMap<XorEncodingSingle>,
}
impl XorEncoding {
	pub fn random_from_tyi<R: Rng>(
		rng: &mut R,
		tyi_output_coding: RowMap<XorEncodingSingle>,
	) -> Self {
		let xor_high: PurposeMap<XorEncodingSingle> = PurposeMap([
			tyi_output_coding[SRow::_0].clone(),
			tyi_output_coding[SRow::_1].clone(),
			rng.random(),
		]);
		let xor_low: PurposeMap<XorEncodingSingle> = PurposeMap([
			tyi_output_coding[SRow::_2].clone(),
			tyi_output_coding[SRow::_3].clone(),
			rng.random(),
		]);
		let xor_output: PurposeMap<XorEncodingSingle> = PurposeMap([
			xor_high[Purpose::Output].clone(),
			xor_low[Purpose::Output].clone(),
			rng.random(),
		]);
		Self {
			xor_high,
			xor_low,
			xor_output,
		}
	}
}

impl XorRound {
	pub fn encode(&mut self, encoding: XorEncoding) -> XorEncodingSingle {
		self.encode_single(Purpose::High, &encoding.xor_high);
		self.encode_single(Purpose::Low, &encoding.xor_low);
		self.encode_single(Purpose::Output, &encoding.xor_output);
		encoding.xor_output[Purpose::Output].clone()
	}
	pub fn encode_single(
		&mut self,
		table_purpose: Purpose,
		encoding: &PurposeMap<XorEncodingSingle>,
	) {
		for i in SColumn::ALL {
			for high_low in HighLow::ALL {
				for j in SColumn::ALL {
					for a in U4::ALL {
						let table = &mut self.0[i][table_purpose][j][high_low][a];
						let output_bijection = NibbleMap(U4::ALL.map(|b| {
							let perm1 = &encoding[Purpose::High].0[i][j][high_low];
							let a = perm1.unmap(a);
							let perm2 = &encoding[Purpose::Low].0[i][j][high_low];
							let b = perm2.unmap(b);

							let out = a ^ b;
							let perm3 = &encoding[Purpose::Output].0[i][j][high_low];
							perm3.map(out)
						}));
						*table = Bijection4::new(output_bijection.0);
					}
				}
			}
		}
	}
}
impl Work {
	pub fn encode(
		&mut self,
		input_encoding: &Option<XorEncodingSingle>,
		output_encoding: &Option<RowMap<XorEncodingSingle>>,
		shift: &ShiftRowsBijection,
	) {
		for pos in SPos::all() {
			let tybox_copy = self.0[pos];
			for (x, sv) in self.0[pos].iter_mut() {
				let shifted_index = shift.map(pos);
				let mut temp = x;

				if let Some(input_encoding) = input_encoding {
					temp = input_encoding.unmap(shifted_index, temp);
				}

				let mut res = tybox_copy[temp];

				if let Some(output_encoding) = output_encoding {
					for high_low in HighLow::ALL {
						for column in SColumn::ALL {
							let v = res.column_nibble(column, high_low);

							let v =
								output_encoding[pos.row()].0[pos.column()][column][high_low].map(v);

							res.set_column_nibble(column, high_low, v)
						}
					}
				}

				*sv = res;
			}
		}
	}
}

impl Tbox {
	pub fn encode(
		&mut self,
		input_encoding: &XorEncodingSingle,
		shift_tables: &ShiftRowsBijection,
	) {
		let copy = self.0;
		for (x, out) in self.0.iter_mut() {
			for pos in SPos::all() {
				let shifted_index = shift_tables.map(pos);

				let temp = input_encoding.unmap(shifted_index, x);
				let res = copy[temp][pos];
				out[pos] = res;
			}
		}
	}
}
