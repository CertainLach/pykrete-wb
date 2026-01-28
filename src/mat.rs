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

use rand::Rng;
use tracing::debug;

use crate::{Column, Row};

fn int2vec(iin: u32) -> VecGF2<32> {
	let mut out = [false; 32];
	for i in 0..32 {
		out[31 - i] = ((iin >> i) & 1) == 1;
	}
	VecGF2::from_array(out)
}
fn vec2int(vin: VecGF2<32>) -> u32 {
	let mut out = 0;
	for i in vin.data.iter() {
		out *= 2;
		if *i == GF2(true) {
			out += 1;
		}
	}
	out
}
fn i82vec(iin: u8) -> VecGF2<8> {
	let mut out = [false; 8];
	for i in 0..8 {
		out[7 - i] = ((iin >> i) & 1) == 1;
	}
	VecGF2::from_array(out)
}
fn vec2i8(vin: VecGF2<8>) -> u8 {
	let mut out = 0;
	for i in vin.data.iter() {
		out *= 2;
		if *i == GF2(true) {
			out += 1;
		}
	}
	out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GF2(pub bool);

impl GF2 {
	pub const ZERO: Self = GF2(false);
	pub const ONE: Self = GF2(true);

	#[inline]
	pub fn new(value: bool) -> Self {
		GF2(value)
	}
}

impl std::ops::Add for GF2 {
	type Output = Self;
	#[inline]
	fn add(self, rhs: Self) -> Self::Output {
		GF2(self.0 ^ rhs.0)
	}
}

impl std::ops::AddAssign for GF2 {
	#[inline]
	fn add_assign(&mut self, rhs: Self) {
		self.0 ^= rhs.0;
	}
}

impl std::ops::Mul for GF2 {
	type Output = Self;
	#[inline]
	fn mul(self, rhs: Self) -> Self::Output {
		GF2(self.0 & rhs.0)
	}
}

impl fmt::Display for GF2 {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", if self.0 { "1" } else { "0" })
	}
}

impl From<bool> for GF2 {
	fn from(b: bool) -> Self {
		GF2(b)
	}
}

#[derive(Clone)]
pub struct MatGF2<const N: usize> {
	data: [[GF2; N]; N],
}

impl MatGF2<8> {
	pub fn identity() -> Self {
		let mut m = Self::new();
		for i in 0..8 {
			m.set(i, i, GF2::ONE);
		}
		m
	}

	pub fn from_column_bytes(cols: [u8; 8]) -> Self {
		let mut m = Self::new();
		for col in 0..8 {
			for row in 0..8 {
				m.set(7 - row, 7 - col, GF2((cols[col] >> row) & 1 == 1));
			}
		}
		m
	}

	pub fn to_column_bytes(&self) -> [u8; 8] {
		let mut cols = [0u8; 8];
		for col in 0..8 {
			for row in 0..8 {
				if self.get(7 - row, 7 - col) == GF2::ONE {
					cols[col] |= 1 << row;
				}
			}
		}
		cols
	}

	pub fn apply(&self, v: u8) -> u8 {
		vec2i8(*self * i82vec(v))
	}

	pub fn apply_inverse(&self, v: u8) -> u8 {
		self.inverse().expect("inversible expected").apply(v)
	}

	pub fn compose(&self, other: &Self) -> Self {
		let mut result = Self::new();
		for i in 0..8 {
			for j in 0..8 {
				let mut sum = GF2::ZERO;
				for k in 0..8 {
					sum = sum + self.get(i, k) * other.get(k, j);
				}
				result.set(i, j, sum);
			}
		}
		result
	}
}
impl MatGF2<32> {
	pub fn apply_row(&self, v: Column) -> Column {
		Column::from_bytes(vec2int(*self * int2vec(u32::from_be_bytes(v.as_bytes()))).to_be_bytes())
	}
	pub fn apply_inverse_row(&self, v: Column) -> Column {
		self.inverse()
			.expect("inversible matrix expected")
			.apply_row(v)
	}
	pub fn apply_column(&self, v: Row) -> Row {
		Row::from_bytes(vec2int(*self * int2vec(u32::from_be_bytes(v.as_bytes()))).to_be_bytes())
	}
	pub fn apply_inverse_column(&self, v: Row) -> Row {
		self.inverse()
			.expect("inversible matrix expected")
			.apply_column(v)
	}
	pub fn concat(mat: [MatGF2<8>; 4]) -> Self {
		let mut concat = Self::new();
		for i in 0..8 {
			for j in 0..8 {
				for (p, off) in mat.iter().zip([0, 8, 16, 24].into_iter()) {
					concat.set(i + off, j + off, p.get(i, j));
				}
			}
		}
		concat
	}
}

impl<const N: usize> MatGF2<N> {
	pub fn random_invertible(rng: &mut impl Rng) -> Self {
		let mut out = MatGF2::random(rng);
		let mut i = 1;
		while out.determinant() == GF2(false) {
			if i % 10 == 0 {
				debug!("matrix #{i} was not invertible");
			}
			out = MatGF2::random(rng);
			i += 1;
		}
		debug!("invertible matrix created in {i} retries");
		out
	}
	pub fn random(rnd: &mut impl Rng) -> Self {
		let mut out = Self::new();
		for ele in out.data.iter_mut() {
			for ele in ele.iter_mut() {
				*ele = GF2(rnd.random())
			}
		}
		out
	}
	pub const fn new() -> Self {
		Self {
			data: [[GF2::ZERO; N]; N],
		}
	}

	#[inline]
	pub fn get(&self, row: usize, col: usize) -> GF2 {
		self.data[row][col]
	}

	#[inline]
	pub fn set(&mut self, row: usize, col: usize, value: GF2) {
		self.data[row][col] = value;
	}

	pub fn determinant(&self) -> GF2 {
		if N == 0 {
			return GF2::ONE;
		}

		let mut matrix = *self;
		let mut det = GF2::ONE;

		for i in 0..N {
			// Find pivot
			let mut pivot_row = None;
			for k in i..N {
				if matrix.data[k][i] == GF2::ONE {
					pivot_row = Some(k);
					break;
				}
			}

			let pivot_row = match pivot_row {
				Some(r) => r,
				None => return GF2::ZERO,
			};

			// Swap rows if needed
			if pivot_row != i {
				matrix.data.swap(i, pivot_row);
				det += GF2::ONE;
			}

			// Eliminate below pivot
			for k in (i + 1)..N {
				if matrix.data[k][i] == GF2::ONE {
					for j in 0..N {
						matrix.data[k][j] += matrix.data[i][j];
					}
				}
			}
		}

		det
	}

	pub fn inverse(&self) -> Option<Self> {
		if N == 0 {
			return Some(*self);
		}

		let mut aug = [[[GF2::ZERO; 2]; N]; N];

		// Create augmented matrix [A | I]
		for i in 0..N {
			for j in 0..N {
				aug[i][j][0] = self.data[i][j]; // Left part: original matrix
				aug[i][j][1] = if i == j { GF2::ONE } else { GF2::ZERO }; // Right part: identity
			}
		}

		// Gaussian elimination
		for i in 0..N {
			// Find pivot
			let mut pivot_row = None;
			for k in i..N {
				if aug[k][i][0] == GF2::ONE {
					pivot_row = Some(k);
					break;
				}
			}

			// Matrix is singular
			let pivot_row = pivot_row?;

			// Swap rows if needed
			if pivot_row != i {
				aug.swap(i, pivot_row);
			}

			// Eliminate column
			for k in 0..N {
				if k != i && aug[k][i][0] == GF2::ONE {
					for j in 0..N {
						aug[k][j][0] += aug[i][j][0];
						aug[k][j][1] += aug[i][j][1];
					}
				}
			}
		}

		// Extract inverse from right half
		let mut result = Self::new();
		for i in 0..N {
			for j in 0..N {
				result.data[i][j] = aug[i][j][1];
			}
		}

		Some(result)
	}
	pub fn from_array(data: [[bool; N]; N]) -> Self {
		let mut result = Self::new();
		for i in 0..N {
			for j in 0..N {
				result.data[i][j] = GF2(data[i][j])
			}
		}
		result
	}
}

impl<const N: usize> Default for MatGF2<N> {
	fn default() -> Self {
		Self::new()
	}
}

impl<const N: usize> Copy for MatGF2<N> {}

impl<const N: usize> fmt::Display for MatGF2<N> {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		for i in 0..N {
			write!(f, "[")?;
			for j in 0..N {
				write!(f, "{}", self.data[i][j])?;
				if j < N - 1 {
					write!(f, " ")?;
				}
			}
			writeln!(f, "]")?;
		}
		Ok(())
	}
}

#[derive(Clone)]
pub struct VecGF2<const N: usize> {
	pub data: [GF2; N],
}

impl<const N: usize> VecGF2<N> {
	pub const fn new() -> Self {
		Self {
			data: [GF2::ZERO; N],
		}
	}

	pub fn from_array(data: [bool; N]) -> Self {
		let mut result = Self::new();
		for i in 0..N {
			result.data[i] = GF2(data[i]);
		}
		result
	}

	#[inline]
	pub fn get(&self, index: usize) -> GF2 {
		self.data[index]
	}

	#[inline]
	pub fn set(&mut self, index: usize, value: GF2) {
		self.data[index] = value;
	}

	pub const fn len(&self) -> usize {
		N
	}
}

impl<const N: usize> Copy for VecGF2<N> {}

impl<const N: usize> fmt::Display for VecGF2<N> {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "[")?;
		for i in 0..N {
			write!(f, "{}", self.data[i])?;
			if i < N - 1 {
				write!(f, " ")?;
			}
		}
		write!(f, "]")
	}
}

impl<const N: usize> std::ops::Mul<VecGF2<N>> for MatGF2<N> {
	type Output = VecGF2<N>;

	fn mul(self, rhs: VecGF2<N>) -> Self::Output {
		let mut result = VecGF2::new();
		for i in 0..N {
			let mut sum = GF2::ZERO;
			for j in 0..N {
				sum += self.data[i][j] * rhs.data[j];
			}
			result.data[i] = sum;
		}
		result
	}
}

#[cfg(test)]
mod tests {
	use rand::rng;

	use super::*;

	#[test]
	fn gf2_math() {
		assert_eq!(GF2::ZERO + GF2::ZERO, GF2::ZERO);
		assert_eq!(GF2::ZERO + GF2::ONE, GF2::ONE);
		assert_eq!(GF2::ONE + GF2::ONE, GF2::ZERO);
		assert_eq!(GF2::ONE * GF2::ONE, GF2::ONE);
		assert_eq!(GF2::ZERO * GF2::ONE, GF2::ZERO);
	}

	#[test]
	fn det_2x2() {
		let matrix = MatGF2::<2>::from_array([[true, false], [false, true]]);
		assert_eq!(matrix.determinant(), GF2::ONE);

		let matrix = MatGF2::<2>::from_array([[true, true], [true, true]]);
		assert_eq!(matrix.determinant(), GF2::ZERO);
	}

	#[test]
	fn det_3x3() {
		let matrix = MatGF2::<3>::from_array([
			[true, false, false],
			[false, true, false],
			[false, false, true],
		]);
		assert_eq!(matrix.determinant(), GF2::ONE);
	}

	#[test]
	fn det_const() {
		const MATRIX: MatGF2<2> = MatGF2::new();
		assert_eq!(MATRIX.determinant(), GF2::ZERO);
	}

	#[test]
	fn random_invertible() {
		let m = <MatGF2<32>>::random_invertible(&mut rng());
		m.inverse().expect("invertible");
	}

	#[test]
	fn mul_by_vec() {
		let matrix = MatGF2::<3>::from_array([
			[true, false, true],
			[false, true, false],
			[true, true, false],
		]);
		let vec = VecGF2::<3>::from_array([true, true, false]);
		let result = matrix * vec;

		// Expected: [1*1+0*1+1*0, 0*1+1*1+0*0, 1*1+1*1+0*0] = [1, 1, 0]
		assert_eq!(result.get(0), GF2::ONE);
		assert_eq!(result.get(1), GF2::ONE);
		assert_eq!(result.get(2), GF2::ZERO);
	}

	#[test]
	fn identity_vec() {
		let identity = MatGF2::<3>::from_array([
			[true, false, false],
			[false, true, false],
			[false, false, true],
		]);
		let vec = VecGF2::<3>::from_array([true, false, true]);
		let result = identity * vec;

		assert_eq!(result.get(0), vec.get(0));
		assert_eq!(result.get(1), vec.get(1));
		assert_eq!(result.get(2), vec.get(2));
	}
}
