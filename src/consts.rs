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

use crate::dual::{Poly};
use crate::xor::{Bijection4, ShiftRowsBijection};
use crate::{State, StateMap, U4};

pub const SHIFT_ROWS_TAB: ShiftRowsBijection = ShiftRowsBijection(Bijection4::new([
	U4::_0,
	U4::_5,
	U4::_10,
	U4::_15,
	U4::_4,
	U4::_9,
	U4::_14,
	U4::_3,
	U4::_8,
	U4::_13,
	U4::_2,
	U4::_7,
	U4::_12,
	U4::_1,
	U4::_6,
	U4::_11,
]));
pub(crate) const INV_SHIFT_ROWS_TAB: ShiftRowsBijection =
	ShiftRowsBijection(SHIFT_ROWS_TAB.0.invert());

// Finite field GF(2^8) multiplication of a and b
pub const fn gf_mul_pol_slow(a: u8, b: u8, pol: Poly) -> u8 {
	let pa = Poly(a as u16);
	let pb = Poly(b as u16);
	let prod = pa.mulmod(pb, pol);
	let prod = prod.0;
	assert!(prod <= 0xff);
	prod as u8
}

// Performs the ShiftRows step. All rows are shifted cylindrically to the left.
pub fn shift_rows(state: &mut State, tab: &ShiftRowsBijection) {
	let copy = *state;

	*state = State(StateMap::from_fn(|pos| copy.0[tab.map(pos)]));
}
