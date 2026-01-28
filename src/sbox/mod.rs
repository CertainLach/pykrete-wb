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

use crate::dual::{AFFINE_CONSTANT, AFFINE_MATRIX, Dual, IRREDUCIBLE_POLYNOMIALS, Poly, Q};
use crate::{X, XArr};

use self::generator::random_sbox;

mod generator;

pub trait SBox {
	fn sub_byte(&self, v: X) -> X;
	fn sub_word(&self, v: u32) -> u32 {
		u32::from_ne_bytes(v.to_ne_bytes().map(|v| self.sub_byte(X(v)).0))
	}
}
impl<S> SBox for &S
where
	S: SBox,
{
	fn sub_byte(&self, v: X) -> X {
		(*self).sub_byte(v)
	}
}

#[derive(Debug)]
pub struct PrecomputedSBox(XArr<X>);
impl PrecomputedSBox {
	pub fn aes_standard() -> Self {
		Self::for_dual(Dual::STANDARD)
	}

	/// Note: current sbox generation algorithm is not great
	/// It rejects obviously bad sboxes, but the set of checked sbox properties is not great
	pub fn random(rng: &mut impl Rng) -> Self {
		Self(random_sbox(rng))
	}
	pub fn inversed(v: &impl SBox) -> Self {
		let mut out = [X(0); 256];

		for i in 0..=255 {
			out[v.sub_byte(X(i)).0 as usize] = X(i);
		}

		Self(XArr(out))
	}
	pub fn inverse(&self) -> Self {
		Self::inversed(self)
	}

	pub fn for_dual(dual: Dual) -> Self {
		fn dual_sbox_byte(x: X, q: Q, q_inv: Q) -> X {
			let y = q_inv.apply(x);

			// Multiplicative inverse in standard (IRREDUCIBLE_POLYNOMIALS[0]) field
			// After Q^{-1}, y is in standard representation
			let inv = if y.0 == 0 {
				0
			} else {
				let p = Poly(y.0 as u16);
				let ip = p.inv(IRREDUCIBLE_POLYNOMIALS[0]);
				ip.0 as u8
			};

			q.apply(X(AFFINE_MATRIX.apply(inv) ^ AFFINE_CONSTANT))
		}

		let q = Q::for_dual(dual);
		let q_inv = q.inverse();
		Self(XArr::from_fn_par(|x| dual_sbox_byte(x, q, q_inv)))
	}
}

impl SBox for PrecomputedSBox {
	fn sub_byte(&self, v: X) -> X {
		self.0[v]
	}
}

pub fn invert_sbox(sbox: &[u8; 256]) -> [u8; 256] {
	let mut inv = [0u8; 256];
	for (i, &val) in sbox.iter().enumerate() {
		inv[val as usize] = i as u8;
	}
	inv
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn standard_sbox_is_standard() {
		let sbox = PrecomputedSBox::aes_standard();

		assert_eq!(sbox.sub_byte(X(0x00)), X(0x63));
		assert_eq!(sbox.sub_byte(X(0x01)), X(0x7c));
	}
}
