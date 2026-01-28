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

use std::fmt;
use std::sync::LazyLock;

use rand::Rng;
use rand::seq::IndexedRandom;

use crate::consts::gf_mul_pol_slow;
use crate::karroumi::Delta;
use crate::key::Key;
use crate::mat::MatGF2;
use crate::tbox::Tbox;
use crate::tybox::Work;
use crate::{SPos, State, Word, X, XArr};

#[derive(Clone, Copy, Debug)]
pub struct Dual {
	pub poly: Poly,
	pub square_power: usize,
}

impl Dual {
	pub const STANDARD: Self = Self {
		poly: IRREDUCIBLE_POLYNOMIALS[0],
		square_power: 0,
	};

	pub const fn new(poly: Poly, square_power: usize) -> Self {
		assert!(square_power <= 7);
		Self { poly, square_power }
	}

	pub fn random<R: Rng>(rng: &mut R) -> Self {
		let poly = IRREDUCIBLE_POLYNOMIALS.choose(rng).expect("not empty");
		let square_power = rng.random_range(0..8);
		Self::new(*poly, square_power)
	}

	pub fn delta(self, to: Dual) -> Delta {
		let from = Q::for_dual(self);
		let to = Q::for_dual(to);
		from.delta(to)
	}
}

#[derive(Clone, Copy)]
pub struct Q(MatGF2<8>);
impl Q {
	pub fn for_dual(dual: Dual) -> Self {
		let field_iso = IRREDUCIBLE_POLYNOMIALS[0].isomorphism(dual.poly);

		if dual.square_power == 0 {
			return Self(MatGF2::identity().compose(&field_iso));
		}
		let q = dual.poly.squaring_matrix();
		let mut result = q;
		let mut i = 1;
		while i < dual.square_power {
			result = result.compose(&q);
			i += 1;
		}
		Self(result.compose(&field_iso))
	}

	pub fn delta(&self, to: Q) -> Delta {
		let q_from_inv = self.0.inverse().expect("Q is invertible");
		let matrix = to.0.compose(&q_from_inv);

		Delta::new(matrix)
	}

	pub fn inverse(&self) -> Q {
		Q(self.0.inverse().expect("Q is inversible"))
	}

	pub fn apply(&self, x: X) -> X {
		X(self.0.apply(x.0))
	}
}

impl State {
	pub fn apply_q(&mut self, q: Q) {
		for pos in SPos::all() {
			self[pos] = q.apply(self[pos]);
		}
	}
	pub fn apply_q_inv(&mut self, q: Q) {
		self.apply_q(q.inverse());
	}
}
impl Work {
	pub fn apply_q(&mut self, q: Q) {
		for pos in crate::SPos::all() {
			let old = self.0[pos];
			self.0[pos] = XArr::from_fn(|x| old[q.apply(x)]);
		}
	}
}
impl Tbox {
	pub fn apply_q_inv(&mut self, q: Q) {
		let q_inv = q.inverse();
		for x in X::all() {
			self.0[x].apply_q(q_inv);
		}
	}
}
impl<const NK: usize> Key<NK> {
	pub fn apply_q(&mut self, q: Q) {
		let mut out = [Word([0, 0, 0, 0]); NK];
		for (i, word) in self.0.iter().enumerate() {
			let transformed = word.to_bytes().map(|b| q.0.apply(b));
			out[i] = Word::from_bytes(transformed);
		}
		self.0 = out;
	}
}

pub fn find_root(source_poly: Poly, target_poly: Poly) -> Option<u8> {
	for alpha in 1u8..=255 {
		let mut val = 0u8;
		let mut alpha_power = 1u8;
		for i in 0..9 {
			if source_poly.has_pow(i) {
				val ^= alpha_power;
			}
			if i < 8 {
				alpha_power = gf_mul_pol_slow(alpha_power, alpha, target_poly);
			}
		}
		if val == 0 {
			return Some(alpha);
		}
	}
	None
}

/// c' = Q^q(phi(c)), where phi is the field isomorphism from IRREDUCIBLE_POLYNOMIAL[0] to poly
/// Ensures transform(c * x) = c' * transform(x)
pub fn transform_coeff_full(coeff: u8, dual: Dual) -> u8 {
	let coeff_in_poly = IRREDUCIBLE_POLYNOMIALS[0]
		.isomorphism(dual.poly)
		.apply(coeff);

	// Frobenius power in poly
	let mut result = coeff_in_poly;
	for _ in 0..dual.square_power {
		result = Poly(result as u16).mulmod(Poly(result as u16), dual.poly).0 as u8;
	}
	result
}

#[derive(Clone, Copy, Debug)]
pub struct MixColCoeffs {
	pub c2: u8,
	pub c3: u8,
	pub c9: u8,
	pub c11: u8,
	pub c13: u8,
	pub c14: u8,
}

impl MixColCoeffs {
	pub const STANDARD: Self = Self {
		c2: 2,
		c3: 3,
		c9: 9,
		c11: 11,
		c13: 13,
		c14: 14,
	};

	pub fn for_dual(dual: Dual) -> Self {
		Self {
			c2: transform_coeff_full(2, dual),
			c3: transform_coeff_full(3, dual),
			c9: transform_coeff_full(9, dual),
			c11: transform_coeff_full(11, dual),
			c13: transform_coeff_full(13, dual),
			c14: transform_coeff_full(14, dual),
		}
	}
}

pub const IRREDUCIBLE_POLYNOMIALS: [Poly; 30] = {
	let mut res = [Poly(0); 30];
	let mut found = 0;

	let mut v = 0x100;
	while v <= 0x1ff {
		let p = Poly(v);
		if p.is_irreducible() {
			res[found] = p;
			found += 1;
		}
		v += 1;
	}

	assert!(found == 30, "there is exactly 30 irreducible polynomials");
	res
};
// Conway polynomial, I don't quite understand the framework here, values taken from conway.pdf
pub const INITIAL_GENERATORS: [u8; 30] = [
	0x03, 0x02, 0x49, 0x2A, 0x2A, 0x4E, 0x0D, 0x0B, 0x2F, 0x2A, 0x21, 0x17, 0x16, 0x26, 0x07, 0x32,
	0x21, 0x06, 0x32, 0x37, 0x1A, 0x1C, 0x25, 0x21, 0x1B, 0x1D, 0x06, 0x21, 0x21, 0x07,
];

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Poly(pub u16);
impl fmt::Debug for Poly {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "Poly(")?;
		let mut had = false;
		for i in (0..16).rev() {
			let has = self.has_pow(i);
			if !has {
				continue;
			}
			if had {
				write!(f, " + ")?;
			}
			had = true;
			if i == 1 {
				write!(f, "x")?;
			} else if i == 0 {
				write!(f, "1")?;
			} else {
				write!(f, "x.pow({i})")?;
			}
		}
		write!(f, ")")
	}
}
impl Poly {
	const fn has_pow(&self, pow: usize) -> bool {
		assert!(pow < 16);
		(self.0 >> pow) & 1 == 1
	}
	const fn toggle_pow(&mut self, pow: usize) {
		assert!(pow < 16);
		self.0 ^= 1 << pow
	}
	pub const fn set_pow(&mut self, pow: usize) {
		assert!(pow < 16);
		self.0 |= 1 << pow
	}
	const fn deg(&self) -> Option<usize> {
		let mut i = 16;
		while i > 0 {
			i -= 1;
			if self.has_pow(i) {
				return Some(i);
			}
		}
		None
	}
	fn squaring_matrix(&self) -> MatGF2<8> {
		let mut cols = [0u8; 8];
		for i in 0..8 {
			let basis = 1u8 << i;
			let squared = Poly(basis as u16).mulmod(Poly(basis as u16), *self);
			cols[i] = squared.0 as u8;
		}
		MatGF2::from_column_bytes(cols)
	}
	fn isomorphism(&self, target: Poly) -> MatGF2<8> {
		if *self == target {
			return MatGF2::identity();
		}

		let root = find_root(*self, target).expect("irreducible polynomials have root");

		// i is the image of the basis element x^i, which is root^i in the target_poly field
		let mut cols = [0u8; 8];
		let mut root_power = 1u8;
		for i in 0..8 {
			cols[i] = root_power;
			if i < 7 {
				root_power = gf_mul_pol_slow(root_power, root, target);
			}
		}
		MatGF2::from_column_bytes(cols)
	}
	pub const fn add(&self, rhs: Self) -> Self {
		Self(self.0 ^ rhs.0)
	}
	pub const fn mul(&self, rhs: Self) -> Self {
		let mut r = Poly(0);
		let mut i = 0;
		while i < 16 {
			if self.has_pow(i) {
				let mut j = 0;
				while j < 16 {
					if rhs.has_pow(j) {
						r.toggle_pow(i + j);
					}
					j += 1;
				}
			}
			i += 1;
		}
		r
	}
	pub const fn modulus(self, modulus: Poly) -> Poly {
		let (_, rem) = self.divmod(modulus);
		rem
	}
	pub const fn mulmod(&self, b: Poly, modulus: Poly) -> Poly {
		self.mul(b).modulus(modulus)
	}
	pub const fn powmod(self, exp: usize, modulus: Poly) -> Poly {
		let mut result = Poly(0b1);
		let mut power = self;
		let mut e = exp;
		while e > 0 {
			if e & 1 == 1 {
				result = result.mulmod(power, modulus);
			}
			power = power.mulmod(power, modulus);
			e >>= 1;
		}
		result
	}
	pub const fn divmod(self, b: Poly) -> (Poly, Poly) {
		let mut r = self;
		let d = b;
		let Some(deg_d) = d.deg() else {
			panic!("division by zero");
		};
		let mut q = Poly(0);
		while let Some(deg_r) = r.deg()
			&& deg_r >= deg_d
		{
			let shift = deg_r - deg_d;
			q.toggle_pow(shift);
			let mut i = 0;
			while i <= deg_d {
				if d.has_pow(i) {
					r.toggle_pow(shift + i);
				}
				i += 1;
			}
		}
		(q, r)
	}
	pub const fn gcd(self, mut b: Poly) -> Poly {
		let mut a = self;
		while b.deg().is_some() {
			let (_, r) = a.divmod(b);
			a = b;
			b = r;
		}
		a
	}
	pub const fn is_irreducible(self) -> bool {
		let poly = self;
		let Some(n) = poly.deg() else {
			return false;
		};
		if n < 1 {
			return false;
		}
		let x = Poly(0b10);
		let xp = {
			let exp = 1 << n;
			x.powmod(exp, poly)
		};
		let diff = xp.add(x);
		let diff_mod = diff.modulus(poly);
		if diff_mod.0 != 0 {
			return false;
		}

		let mut d = 1;
		while d * d <= n {
			if n.is_multiple_of(d) {
				let k = d;
				{
					// IDENTICAL BLOCKS
					if k < n && k > 0 {
						let xp_k = x.powmod(1 << k, poly);
						let g = xp_k.add(x).gcd(poly);
						if g.0 != 1 {
							return false;
						}
					}
				}
				let k = n / d;
				{
					// IDENTICAL BLOCKS
					if k < n && k > 0 {
						let xp_k = x.powmod(1 << k, poly);
						let g = xp_k.add(x).gcd(poly);
						if g.0 != 1 {
							return false;
						}
					}
				}
			}
			d += 1;
		}
		true
	}
	pub const fn inv(self, modulus: Poly) -> Poly {
		let mut r0 = self;
		let mut r1 = modulus;
		let mut s0: Poly = Poly(1);
		let mut s1: Poly = Poly(0);
		while r1.deg().is_some() {
			let (q, r2) = r0.divmod(r1);
			let s2 = s0.add(q.mul(s1));
			r0 = r1;
			r1 = r2;
			s0 = s1;
			s1 = s2;
		}
		if r0.deg().expect("exists") != 0 || !r0.has_pow(0) {
			panic!("no inv found")
		}
		let mut inv = s0;
		inv = inv.modulus(modulus);
		inv
	}
}

pub struct G([u8; 256], Poly);
impl G {
	pub fn mul(&self, a: X, idx: usize) -> X {
		X(gf_mul_pol_slow(a.0, self.0[idx], self.1))
	}
}

pub const fn make_g(generator: u8, pol: Poly) -> G {
	let mut out = [0; 256];
	let mut i = 0;
	let mut cur = 1;
	while i <= 255 {
		out[i] = cur;
		i += 1;
		cur = gf_mul_pol_slow(cur, generator, pol)
	}
	G(out, pol)
}

pub fn is_generator(alpha: u8, field: Poly) -> bool {
	if !Poly(alpha as u16).is_irreducible() {
		return false;
	}
	for d in [1, 3, 5, 15, 17, 51, 85] {
		if Poly(alpha as u16).powmod(d, field) == Poly(1) {
			return false;
		}
	}
	true
}

pub fn generators_for_pol(pol: Poly, first: u8) -> [u8; 8] {
	let first = first as u16;
	[
		first as u8,
		Poly(first).powmod(2, pol).0 as u8,
		Poly(first).powmod(4, pol).0 as u8,
		Poly(first).powmod(8, pol).0 as u8,
		Poly(first).powmod(16, pol).0 as u8,
		Poly(first).powmod(32, pol).0 as u8,
		Poly(first).powmod(64, pol).0 as u8,
		Poly(first).powmod(128, pol).0 as u8,
	]
}

pub static AFFINE_MATRIX: LazyLock<MatGF2<8>> = LazyLock::new(|| {
	MatGF2::from_column_bytes([
		0b00011111, 0b00111110, 0b01111100, 0b11111000, 0b11110001, 0b11100011, 0b11000111,
		0b10001111,
	])
});

pub const AFFINE_CONSTANT: u8 = 0x63;

pub fn compute_rcon(round: usize, config: Dual) -> u8 {
	if round == 0 {
		return 0;
	}
	let poly = config.poly;
	let base = transform_coeff_full(0x02, config);
	Poly(base as u16).powmod(round - 1, poly).0 as u8
}

#[cfg(test)]
mod tests {
	use crate::dual::{Dual, IRREDUCIBLE_POLYNOMIALS, Q};
	use crate::key::{Aes128Key, RoundKeys, expand_nonstandard_keys};
	use crate::sbox::{PrecomputedSBox, SBox as _};
	use crate::{Security, State, Tables, decrypt_nonstandard, encrypt_nonstandard};
	use rand::rng;
	use test_case::test_case;

	use super::*;

	#[test]
	fn dual() {
		let mut cur = 1;

		for i in 0..255 {
			println!("g[{i}] = {cur:x?}");
			cur = gf_mul_pol_slow(cur, 3, IRREDUCIBLE_POLYNOMIALS[0])
		}
	}

	#[test]
	fn poly_basic_math() {
		assert_eq!(Poly(0b11).mul(Poly(0b101)), Poly(0b1111));
		assert_eq!(Poly(0b101).mul(Poly(0b1010)), Poly(0b100010));

		assert_eq!(Poly(0b101).divmod(Poly(0b11)), (Poly(0b11), Poly(0)));
		assert_eq!(Poly(0b1010).divmod(Poly(0b101)), (Poly(0b10), Poly(0)));
	}

	#[test]
	fn generators() {
		let g = make_g(INITIAL_GENERATORS[0], IRREDUCIBLE_POLYNOMIALS[0]);
		for i in 0..256 {
			println!("g[{i}] = {}", g.0[i]);
		}
	}

	#[test]
	fn squaring_matrix_standard() {
		let poly = IRREDUCIBLE_POLYNOMIALS[0];
		let q = poly.squaring_matrix();

		for x in 0..=255u8 {
			let squared = Poly(x as u16).mulmod(Poly(x as u16), poly);
			let via_matrix = q.apply(x);
			assert_eq!(squared, Poly(via_matrix as u16));
		}
	}

	#[test]
	fn matrix_inverse() {
		let q = IRREDUCIBLE_POLYNOMIALS[0].squaring_matrix();
		let q_inv = q.inverse().expect("Q should be invertible");

		for x in 0..=255u8 {
			let qx = q.apply(x);
			let back = q_inv.apply(qx);
			assert_eq!(back, x, "Q^{{-1}}(Q(x)) == x");
		}
	}

	#[test]
	fn dual_cipher_equivalence() {
		let poly = IRREDUCIBLE_POLYNOMIALS[0];

		let q = poly.squaring_matrix();
		let _q_inv = q.inverse().unwrap();

		let config = Dual::new(poly, 1);
		let dual_sbox = PrecomputedSBox::for_dual(config);
		let standard_sbox = PrecomputedSBox::for_dual(Dual::STANDARD);

		for x in X::all() {
			let standard_out = standard_sbox.sub_byte(x);
			let squared_standard_out = q.apply(standard_out.0);

			let squared_x = q.apply(x.0);
			let dual_out = dual_sbox.sub_byte(X(squared_x));

			assert_eq!(squared_standard_out, dual_out.0, "S^2(Q*x) == Q*S(x)");
		}
	}

	#[test]
	fn mixcol_coeffs_squared() {
		let poly = IRREDUCIBLE_POLYNOMIALS[0];

		let c2_squared = Poly(0x02).mulmod(Poly(0x02), poly);
		assert_eq!(c2_squared, Poly(0x04));

		let c3_squared = Poly(0x03).mulmod(Poly(0x03), poly);
		assert_eq!(c3_squared, Poly(0x05));

		let config = Dual::new(poly, 1);
		let coeffs = MixColCoeffs::for_dual(config);
		assert_eq!(coeffs.c2, 0x04);
		assert_eq!(coeffs.c3, 0x05);
	}

	#[test]
	fn standard_rcon() {
		let expected_rcons: [u8; 10] = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36];
		for (i, &expected) in expected_rcons.iter().enumerate() {
			let computed = compute_rcon(i + 1, Dual::STANDARD);
			assert_eq!(computed, expected);
		}
	}

	#[test]
	fn q8_is_identity() {
		let q = IRREDUCIBLE_POLYNOMIALS[0].squaring_matrix();
		let mut q_power = MatGF2::identity();
		for _ in 0..8 {
			q_power = q_power.compose(&q);
		}

		for x in 0..=255u8 {
			assert_eq!(q_power.apply(x), x, "Frobenius identity");
		}
	}

	#[test_case(false, 0, 0; "standard pol, power")]
	#[test_case(false, 0, 1; "standard pol, non-standard power")]
	#[test_case(false, 1, 0; "non-standard pol, standard power")]
	#[test_case(false, 1, 1; "non-standard pol, power")]
	#[test_case(true, 0, 0; "whitebox, standard pol, power")]
	#[test_case(true, 0, 1; "whitebox, standard pol, non-standard power")]
	#[test_case(true, 1, 0; "whitebox, non-standard pol, standard power")]
	#[test_case(true, 1, 1; "whitebox, non-standard pol, power")]
	fn dual_cipher_equivalence_square(tables: bool, pol: usize, power: usize) {
		let dual = Dual::new(IRREDUCIBLE_POLYNOMIALS[pol], power);

		let key = Aes128Key::KUNG_FU_TEST_VECTOR;
		let standard_rk = key.expand();

		let dual_sbox = PrecomputedSBox::for_dual(dual);
		let dual_rk =
			expand_nonstandard_keys::<_, 4, 9>(&dual_sbox, Aes128Key::KUNG_FU_TEST_VECTOR, dual);

		let mut standard_data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		encrypt_nonstandard(
			&PrecomputedSBox::aes_standard(),
			&mut standard_data,
			&standard_rk,
			Dual::STANDARD,
		);

		let mut dual_data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;

		if tables {
			let mut tables = Tables::from_nonstandard_round_keys(&dual_sbox, &dual_rk, false, dual);
			tables.apply_security(Security::full(), &mut rng());
			tables.cipher(&mut dual_data);
		} else {
			encrypt_nonstandard(&dual_sbox, &mut dual_data, &dual_rk, dual);
		};

		let q = Q::for_dual(dual);
		standard_data.apply_q(q);

		assert_eq!(standard_data, dual_data);

		if tables {
			let mut tables = Tables::from_nonstandard_round_keys(&dual_sbox, &dual_rk, true, dual);
			tables.apply_security(Security::full(), &mut rng());
			tables.cipher(&mut dual_data);
		} else {
			decrypt_nonstandard(&dual_sbox, &mut dual_data, &dual_rk, dual);
		}

		assert_eq!(State::TWO_ONE_NINE_TWO_TEST_VECTOR, dual_data)
	}

	#[test]
	fn all_square_duals() {
		for square_power in 1..8 {
			let dual = Dual::new(IRREDUCIBLE_POLYNOMIALS[0], square_power);

			let key = Aes128Key::new([0x2b; 16]);
			let standard_rk = key.expand();
			let dual_sbox = PrecomputedSBox::for_dual(dual);
			let dual_rk =
				expand_nonstandard_keys::<_, 4, 9>(&dual_sbox, Aes128Key::new([0x2b; 16]), dual);

			let mut standard_data = State::from_bytes([0x32; 16]);
			encrypt_nonstandard(
				&PrecomputedSBox::aes_standard(),
				&mut standard_data,
				&standard_rk,
				Dual::STANDARD,
			);

			// Q*E(P,K)
			let q = Q::for_dual(dual);
			standard_data.apply_q(q);

			let mut dual_data = State::from_bytes([0x32; 16]);
			encrypt_nonstandard(&dual_sbox, &mut dual_data, &dual_rk, dual);

			assert_eq!(standard_data, dual_data, "equivalence fail");
		}
	}

	#[test]
	fn q8_cycles_back() {
		let q = IRREDUCIBLE_POLYNOMIALS[0].squaring_matrix();
		let mut transform = crate::mat::MatGF2::identity();
		for _ in 0..8 {
			transform = transform.compose(&q);
		}

		for x in 0..=255u8 {
			assert_eq!(transform.apply(x), x, "should be identity");
		}
	}

	#[test]
	fn dual_roundtrip_nonwhitebox() {
		let dual = Dual::new(IRREDUCIBLE_POLYNOMIALS[0], 1);

		let dual_rk = expand_nonstandard_keys::<_, 4, 9>(
			&PrecomputedSBox::for_dual(dual),
			Aes128Key::KUNG_FU_TEST_VECTOR,
			dual,
		);
		let dual_sbox = PrecomputedSBox::for_dual(dual);

		let original = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		let mut data = original;

		encrypt_nonstandard(&dual_sbox, &mut data, &dual_rk, dual);
		{
			// Data is different from AES because of Q, apply Q^{-1}
			let q = Q::for_dual(dual);
			let mut standard_data = data.clone();
			standard_data.apply_q_inv(q);

			let mut standard_data_comparison = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
			// Compare with AES
			let normal_rk: RoundKeys<9> = expand_nonstandard_keys(
				&PrecomputedSBox::aes_standard(),
				Aes128Key::KUNG_FU_TEST_VECTOR,
				Dual::STANDARD,
			);
			encrypt_nonstandard(
				&PrecomputedSBox::aes_standard(),
				&mut standard_data_comparison,
				&normal_rk,
				Dual::STANDARD,
			);

			assert_eq!(standard_data_comparison, standard_data);
		}
		decrypt_nonstandard(&dual_sbox, &mut data, &dual_rk, dual);

		assert_eq!(data, original);
	}
}
