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

use crate::tbox::Tboxes;
use crate::ty::Ty;
use crate::{RIArr, Row, SColumn, StateMap, XArr};

#[derive(Default, Clone, Copy, Debug)]
pub struct Work(pub(crate) StateMap<XArr<Row>>);

#[derive(Debug)]
pub struct WorkRounds<const NRM1: usize>(pub(crate) RIArr<Work, NRM1>);
impl<const NRM1: usize> WorkRounds<NRM1> {
	pub fn new_tyi(tboxes: &Tboxes<NRM1>, ty: &Ty) -> Self {
		Self(RIArr::from_fn(|r| {
			Work(StateMap::from_fn(|pos| {
				XArr::from_fn(|x| {
					let i = pos.row();
					let tboxv = tboxes.get(r, x, pos);
					ty.get_row(tboxv, i)
				})
			}))
		}))
	}
	pub fn new_mbl() -> Self {
		let round = Work(StateMap::from_fn(|pos| {
			XArr::from_fn(|x| {
				let mut mblv = Row::default();
				mblv[SColumn(pos.row().0)] = x;
				mblv
			})
		}));
		Self(RIArr::from_fn(|_| round))
	}
}
