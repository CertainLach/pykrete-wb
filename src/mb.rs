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

use crate::{
	Column, ColumnMap, RowMap, SColumn, SPos, SRow, X, XArr, mat::MatGF2, tbox::Tbox, tybox::Work,
	xor::ShiftRowsBijection,
};

pub struct L(RowMap<ColumnMap<MatGF2<8>>>);
impl L {
	pub fn random<R: Rng>(rng: &mut R) -> Self {
		Self(RowMap::from_fn(|_| {
			ColumnMap::from_fn(|_| MatGF2::random_invertible(rng))
		}))
	}
}

impl Work {
	pub fn apply_l_inv(&mut self, l: &L, shift_rows: &ShiftRowsBijection) {
		for pos in SPos::all() {
			let old_tyboxes = self.0[pos];
			for x in X::all() {
				let lv = {
					let pos = shift_rows.map(pos);
					l.0[pos.row()][pos.column()]
				};
				let idx = X(lv.apply_inverse(x.0));
				let v = old_tyboxes[idx];
				self.0[pos][x] = v;
			}
		}
	}
	pub fn apply_l(&mut self, l: &L) {
		for row in SRow::all() {
			for x in X::all() {
				let comp_elem = |out: Column| {
					let concat = MatGF2::<32>::concat(l.0[row].0);
					concat.apply_column(out)
				};

				for column in SColumn::ALL {
					let pos = SPos::row_column(row, column);
					self.0[pos][x] = comp_elem(self.0[pos][x]);
				}
			}
		}
	}
}

impl Tbox {
	pub fn apply_l_inv(&mut self, l: &L, shift_rows: &ShiftRowsBijection) {
		for pos in SPos::all() {
			let old_tboxes_last = XArr::from_fn(|x| self.0[x][pos]);
			for x in X::all() {
				let idx = {
					let pos = shift_rows.map(pos);
					X(l.0[pos.row()][pos.column()].apply_inverse(x.0))
				};
				self.0[x][pos] = old_tboxes_last[idx];
			}
		}
	}
}

pub struct MB(pub RowMap<MatGF2<32>>);
impl MB {
	pub fn random<R: Rng>(rng: &mut R) -> Self {
		Self(RowMap::from_fn(|_| <MatGF2<32>>::random_invertible(rng)))
	}
}

impl Work {
	pub fn apply_mb(&mut self, mb: &MB) {
		for (ij, xbox) in self.0.iter_mut() {
			let mat = &mb.0[ij.row()];
			for (_, col) in xbox.iter_mut() {
				*col = mat.apply_column(*col)
			}
		}
	}
	pub fn apply_mb_inv(&mut self, mb: &MB) {
		for (ij, xbox) in self.0.iter_mut() {
			let mat = &mb.0[ij.row()];
			for (_, col) in xbox.iter_mut() {
				*col = mat.apply_inverse_column(*col)
			}
		}
	}
}
