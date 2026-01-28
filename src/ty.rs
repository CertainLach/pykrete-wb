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

use crate::consts::gf_mul_pol_slow;
use crate::dual::{Dual, MixColCoeffs, Poly};
use crate::{Column, ColumnMap, SColumn, X, XArr};

pub struct Ty(ColumnMap<XArr<ColumnMap<X>>>);
impl Ty {
	pub fn new(inv: bool, config: Dual) -> Self {
		let coeffs = MixColCoeffs::for_dual(config);
		let poly = config.poly;

		if inv {
			Self::new_inv(&coeffs, poly)
		} else {
			Self::new_forward(&coeffs, poly)
		}
	}

	fn gf_mul(x: X, coeff: u8, poly: Poly) -> X {
		X(gf_mul_pol_slow(x.0, coeff, poly))
	}

	fn new_forward(coeffs: &MixColCoeffs, poly: Poly) -> Self {
		let mut out: ColumnMap<XArr<ColumnMap<X>>> = Default::default();
		for x in X::all() {
			let m2 = Self::gf_mul(x, coeffs.c2, poly);
			let m3 = Self::gf_mul(x, coeffs.c3, poly);

			out[SColumn::_0][x][SColumn::_0] = m2;
			out[SColumn::_1][x][SColumn::_0] = m3;
			out[SColumn::_2][x][SColumn::_0] = x;
			out[SColumn::_3][x][SColumn::_0] = x;

			out[SColumn::_0][x][SColumn::_1] = x;
			out[SColumn::_1][x][SColumn::_1] = m2;
			out[SColumn::_2][x][SColumn::_1] = m3;
			out[SColumn::_3][x][SColumn::_1] = x;

			out[SColumn::_0][x][SColumn::_2] = x;
			out[SColumn::_1][x][SColumn::_2] = x;
			out[SColumn::_2][x][SColumn::_2] = m2;
			out[SColumn::_3][x][SColumn::_2] = m3;

			out[SColumn::_0][x][SColumn::_3] = m3;
			out[SColumn::_1][x][SColumn::_3] = x;
			out[SColumn::_2][x][SColumn::_3] = x;
			out[SColumn::_3][x][SColumn::_3] = m2;
		}
		Self(out)
	}

	fn new_inv(coeffs: &MixColCoeffs, poly: Poly) -> Self {
		let mut out: ColumnMap<XArr<ColumnMap<X>>> = Default::default();
		for x in X::all() {
			let m9 = Self::gf_mul(x, coeffs.c9, poly);
			let m11 = Self::gf_mul(x, coeffs.c11, poly);
			let m13 = Self::gf_mul(x, coeffs.c13, poly);
			let m14 = Self::gf_mul(x, coeffs.c14, poly);

			out[SColumn::_0][x][SColumn::_0] = m14;
			out[SColumn::_1][x][SColumn::_0] = m11;
			out[SColumn::_2][x][SColumn::_0] = m13;
			out[SColumn::_3][x][SColumn::_0] = m9;

			out[SColumn::_0][x][SColumn::_1] = m9;
			out[SColumn::_1][x][SColumn::_1] = m14;
			out[SColumn::_2][x][SColumn::_1] = m11;
			out[SColumn::_3][x][SColumn::_1] = m13;

			out[SColumn::_0][x][SColumn::_2] = m13;
			out[SColumn::_1][x][SColumn::_2] = m9;
			out[SColumn::_2][x][SColumn::_2] = m14;
			out[SColumn::_3][x][SColumn::_2] = m11;

			out[SColumn::_0][x][SColumn::_3] = m11;
			out[SColumn::_1][x][SColumn::_3] = m13;
			out[SColumn::_2][x][SColumn::_3] = m9;
			out[SColumn::_3][x][SColumn::_3] = m14;
		}
		Self(out)
	}

	pub fn get_column(&self, tboxv: X, i: SColumn) -> Column {
		Column(self.0[i][tboxv].0)
	}
}
