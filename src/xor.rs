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
use rand::distr::StandardUniform;
use rand::prelude::*;

use crate::ColumnMap;
use crate::HighLow;
use crate::HighLowMap;
use crate::Purpose;
use crate::PurposeMap;
use crate::SColumn;
use crate::SPos;
use crate::{NibbleMap, RI, RIArr, U4};

#[derive(Clone, Copy)]
pub struct ShiftRowsBijection(pub(crate) Bijection4);
impl ShiftRowsBijection {
	pub(crate) const IDENTITY: Self = Self(Bijection4::IDENTITY);
	pub fn map(&self, pos: SPos) -> SPos {
		SPos(self.0.map(pos.0))
	}
}

#[derive(Clone, Copy)]
pub struct Bijection4([u8; 8]);

impl Distribution<Bijection4> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Bijection4 {
		let mut identity = U4::ALL;
		identity.shuffle(rng);
		Bijection4::new(identity)
	}
}

impl Bijection4 {
	pub(crate) const IDENTITY: Self = Self::new(U4::ALL);
	fn xor_identity(i: U4) -> Self {
		Bijection4::new(U4::ALL.map(|j| i ^ j))
	}
	pub(crate) const fn new(v: [U4; 16]) -> Self {
		let mut out = [0u8; 8];
		let mut i = 0;
		let mut encountered = 0u16;
		while i < 8 {
			let a = v[i * 2].as_index();
			let b = v[i * 2 + 1].as_index();
			encountered |= 1 << a;
			encountered |= 1 << b;

			out[i] = ((a << 4) | b) as u8;
			i += 1;
		}
		assert!(encountered.count_ones() == 16);
		Self(out)
	}
	pub const fn map(&self, i: U4) -> U4 {
		let i = i.as_index();
		let ii = i >> 1;
		let b = self.0[ii];
		U4::from_index(if i & 1 == 1 { b & 0xf } else { b >> 4 } as usize)
	}
	pub fn unmap(&self, v: U4) -> U4 {
		for i in U4::ALL {
			if self.map(i) == v {
				return i;
			}
		}
		unreachable!()
	}
	pub fn mapunmap(&self, v: U4, inv: bool) -> U4 {
		if inv { self.unmap(v) } else { self.map(v) }
	}
	pub const fn invert(self) -> Self {
		let mut out = [U4::_0; 16];

		let mut i = 0;
		while i < 16 {
			out[self.map(U4::from_index(i)).as_index()] = U4::from_index(i);
			i += 1;
		}

		Self::new(out)
	}
}
impl fmt::Debug for Bijection4 {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		let v = U4::ALL.map(|i| self.map(i));
		write!(f, "Bijection(")?;
		for i in v {
			write!(f, "{i:?}")?;
		}
		write!(f, ")")
	}
}
impl Default for Bijection4 {
	fn default() -> Self {
		Self::IDENTITY
	}
}

pub struct XorPartial<'t>(&'t PurposeMap<ColumnMap<HighLowMap<NibbleMap<Bijection4>>>>);
impl<'t> XorPartial<'t> {
	pub(crate) fn map(
		&self,
		purpose: Purpose,
		column: SColumn,
		high_low: HighLow,
		a: U4,
		b: U4,
	) -> U4 {
		self.0[purpose][column][high_low][a].map(b)
	}
}

#[derive(Debug, Clone, Copy)]
pub struct XorRound(pub(crate) ColumnMap<PurposeMap<ColumnMap<HighLowMap<NibbleMap<Bijection4>>>>>);
impl XorRound {
	fn identity() -> Self {
		let nib = NibbleMap::from_fn(Bijection4::xor_identity);
		Self(ColumnMap::from_fn(|_| {
			PurposeMap::from_fn(|_| ColumnMap::from_fn(|_| HighLowMap::from_fn(|_| nib)))
		}))
	}
}

#[derive(Debug)]
pub struct Xor<const NRM1: usize>(pub RIArr<XorRound, NRM1>);

impl<const NRM1: usize> Xor<NRM1> {
	pub(crate) fn partial_map(&self, r: RI, j: SColumn) -> XorPartial<'_> {
		XorPartial(&self.0[r].0[j])
	}

	pub fn identity() -> Self {
		let round = XorRound::identity();
		Self(RIArr::from_fn(|_| round))
	}
}
