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
use std::mem::{MaybeUninit, transmute};
use std::ops::{BitXor, BitXorAssign, Index, IndexMut};
use std::{array, fmt};

use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use zeroize::Zeroize;

use crate::dual::{Dual, MixColCoeffs, Q};
use crate::internal::XorEncoding;
use crate::mb::{L, MB};

use self::consts::{INV_SHIFT_ROWS_TAB, SHIFT_ROWS_TAB, shift_rows};

use self::internal::XorEncodingSingle;
use self::key::{Aes128Key, Aes192Key, Aes256Key, Key, RoundKeys};
use self::tbox::{Tbox, Tboxes};
use self::ty::Ty;
use self::tybox::WorkRounds;
use self::xor::{ShiftRowsBijection, Xor};
mod consts;
pub mod key;
pub mod mat;

pub mod dual;
pub mod encoding;
pub mod hardware;
pub mod internal;
pub mod karroumi;
pub mod mb;
pub mod sbox;
pub mod tbox;
pub mod ty;
pub mod tybox;
pub mod vm;
pub mod xor;

use sbox::{PrecomputedSBox, SBox};

#[derive(Debug, Clone, Copy)]
pub struct Security {
	// Mixing bijections
	// Only enabling one of them is not enough, because it is then pretty easy to remove, but enabling both makes wonders
	// 32-bit bijections
	pub mb: bool,
	// 8-bit bijections
	pub l: bool,
	// Always use mb + l matrices, even if without bijections.
	// Only makes effect with both mb and l are unused.
	pub force_mbl: bool,
	// Always use xor table for mbl.
	// either due to mb/l enabled, or when force_mbl is enabled.
	pub force_mbl_xor: bool,
	/// Always use xor table for tyboxes
	pub force_xor: bool,

	// Internal encodings intermixes values using xor tables, tying rounds together
	pub internal_encodings: bool,
}
impl Security {
	pub fn full() -> Self {
		Self {
			mb: true,
			l: true,
			internal_encodings: true,
			force_mbl: false,
			force_xor: false,
			force_mbl_xor: false,
		}
	}
}

impl<const NRM1: usize> Tables<NRM1> {
	pub fn apply_security<R: Rng>(&mut self, security: Security, rng: &mut R) {
		let Security {
			mb,
			l,
			force_mbl,
			force_xor,
			internal_encodings,
			force_mbl_xor,
		} = security;

		let shift_rows = if self.inv {
			&INV_SHIFT_ROWS_TAB
		} else {
			&SHIFT_ROWS_TAB
		};

		let mbl = &mut self.mbl;
		let mut xor = &mut self.xor;
		let mut xor_mbl = &mut self.xor_mbl;
		let tyboxes = &mut self.tyboxes;
		let tboxes_last = &mut self.tboxes_last;

		if mb {
			for r in RI::all::<NRM1>() {
				let mb = MB::random(rng);
				tyboxes.0[r].apply_mb(&mb);
				mbl.0[r].apply_mb_inv(&mb);
			}
		}
		if l {
			let mut prev_l = None::<L>;

			for r in RI::all::<NRM1>() {
				if let Some(prev_l) = prev_l {
					tyboxes.0[r].apply_l_inv(&prev_l, shift_rows);
				}

				let l = L::random(rng);
				mbl.0[r].apply_l(&l);
				prev_l = Some(l);
			}

			let l = prev_l.expect("at least one round should be processed");
			tboxes_last.apply_l_inv(&l, shift_rows);
		}

		let uses_mbl = mb || l || force_mbl || self.uses_mbl;
		self.uses_mbl |= uses_mbl;

		if internal_encodings {
			let mut prev_round_input_encoding = None::<XorEncodingSingle>;

			for r in RI::all::<NRM1>() {
				for step in [Step::Tybox, Step::Mbl] {
					if matches!(step, Step::Mbl) && !uses_mbl {
						continue;
					}
					let (work, xor, shift) = match step {
						Step::Tybox => {
							self.uses_xor = true;
							(&mut tyboxes.0[r], &mut xor, shift_rows)
						}
						Step::Mbl => {
							self.uses_xor_mbl = true;
							(&mut mbl.0[r], &mut xor_mbl, &ShiftRowsBijection::IDENTITY)
						}
					};
					let tyi_output_coding: RowMap<XorEncodingSingle> = rng.random();

					work.encode(
						&prev_round_input_encoding,
						&Some(tyi_output_coding.clone()),
						shift,
					);

					let xor_encoding = XorEncoding::random_from_tyi(rng, tyi_output_coding);

					prev_round_input_encoding = Some(xor.0[r].encode(xor_encoding));
				}
			}

			if let Some(prev_round_input_encoding) = prev_round_input_encoding {
				tboxes_last.encode(&prev_round_input_encoding, shift_rows);
			}
		}

		if force_mbl_xor {
			self.uses_xor_mbl = true;
		}
		if force_xor {
			self.uses_xor = true;
		}
	}
}

#[derive(Debug)]
pub struct Tables<const NRM1: usize> {
	/// Is that a decryption tables, and inverse shift rows should be used in `cipher`
	pub(crate) inv: bool,

	/// Ty performs MixColumns transform, tybox is Ty+Tbox for every round except last
	///
	/// Consists of
	/// - L^{-1} (optional, 8-bit bijections)
	/// - T
	/// - Ty
	/// - MB (optional, 32-bit bijections)
	///
	/// External input encoding might be applied here for the first round
	pub(crate) tyboxes: WorkRounds<NRM1>,
	/// Are internal encodings used for tyboxes?
	pub(crate) uses_xor: bool,
	/// XOR tables for tyboxes
	pub(crate) xor: Xor<NRM1>,

	/// Inverse of the optional mixing bijections applied to tyboxes for the all rounds except last
	///
	/// Consists of
	/// - MB^{-1} (32-bit bijections)
	/// - L (8-bit bijections)
	///
	/// Application of this table is only required if `uses_mbl` is enabled
	pub(crate) mbl: WorkRounds<NRM1>,
	/// Are MB/L were applied to tyboxes, and it is required to apply MBL lookup table.
	/// Note that the table is always available and populated with identity mixing bijections,
	/// so it should be safe to enable this flag unconditionally.
	pub(crate) uses_mbl: bool,
	/// Are internal encodings used for mbl?
	pub(crate) uses_xor_mbl: bool,
	/// XOR tables for mixed bijections
	pub(crate) xor_mbl: Xor<NRM1>,

	/// Tbox performs AddRoundKeys, SubBytes transform, last round doesn't need anything else
	///
	/// External output encoding might be applied here for the last round
	pub(crate) tboxes_last: Tbox,
}

enum Step {
	Tybox,
	Mbl,
}
impl Step {
	const ALL: [Self; 2] = [Self::Tybox, Self::Mbl];
}

pub struct IV(StateMap<X>);
impl IV {
	#[cfg(test)]
	pub const NIST_CFB_E1: Self = Self::new([
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
		0x0f,
	]);

	pub const fn new(v: [u8; 16]) -> Self {
		Self(unsafe { transmute::<[u8; 16], StateMap<X>>(v) })
	}
}

impl<const NRM1: usize> Tables<NRM1> {
	pub(crate) fn new_base(tyboxes: WorkRounds<NRM1>, tboxes_last: Tbox, inv: bool) -> Self {
		Self {
			tyboxes,
			tboxes_last,
			inv,
			mbl: WorkRounds::new_mbl(),
			uses_mbl: false,
			xor: Xor::identity(),
			xor_mbl: Xor::identity(),
			uses_xor_mbl: false,
			uses_xor: false,
		}
	}

	pub fn from_nonstandard_key<S: SBox + Copy, const NK: usize>(
		sbox: S,
		key: Key<NK>,
		inv: bool,
		dual: Dual,
	) -> Self {
		use self::key::expand_nonstandard_keys;

		Self::from_nonstandard_round_keys(
			sbox,
			&expand_nonstandard_keys(sbox, key, dual),
			inv,
			dual,
		)
	}
	pub fn from_nonstandard_round_keys<S: SBox + Copy>(
		sbox: S,
		round_keys: &RoundKeys<NRM1>,
		inv: bool,
		dual: Dual,
	) -> Self {
		let tboxes = Tboxes::from_round_keys(sbox, round_keys, inv);

		let ty = Ty::new(inv, dual);

		let mut tyboxes = WorkRounds::new_tyi(&tboxes, &ty);

		let mut tboxes_last = tboxes.last;

		let q = Q::for_dual(dual);
		if inv {
			tboxes_last.apply_q_inv(q);
		} else {
			tyboxes.0[RI(0)].apply_q(q);
		}

		Self::new_base(tyboxes, tboxes_last, inv)
	}

	pub fn encrypt_cfb(&self, iv: IV, m: &mut [u8]) {
		let mut cfb_blk = State(iv.0);

		for (i, v) in m.iter_mut().enumerate() {
			if i & 0xf == 0 {
				self.cipher(&mut cfb_blk);
			}
			let pos = SPos::from_index(i & 0xf);
			cfb_blk[pos] ^= X(*v);
			*v = cfb_blk[pos].0;
		}
	}
	pub fn decrypt_cfb(&self, iv: IV, m: &mut [u8]) {
		let mut cfb_blk = State(iv.0);

		for (i, v) in m.iter_mut().enumerate() {
			if i & 0xf == 0 {
				self.cipher(&mut cfb_blk)
			}
			let pos = SPos::from_index(i & 0xf);
			let c = X(*v);
			*v = (cfb_blk[pos] ^ c).0;
			cfb_blk[pos] = c;
		}
	}

	fn cipher(&self, data: &mut State) {
		for r in RI::all::<NRM1>() {
			shift_rows(
				data,
				if self.inv {
					&INV_SHIFT_ROWS_TAB
				} else {
					&SHIFT_ROWS_TAB
				},
			);

			// tbox + ty(i)
			for col in SColumn::all() {
				for step in Step::ALL {
					if matches!(step, Step::Mbl) && !self.uses_mbl {
						continue;
					}
					let work = match step {
						Step::Tybox => &self.tyboxes.0[r],
						Step::Mbl => &self.mbl.0[r],
					};
					let [aa, bb, cc, dd] = SRow::ALL.map(|row| {
						let pos = SPos::column_row(col, row);
						work.0[pos][data[pos]]
					});
					let xor = match step {
						Step::Mbl if self.uses_xor_mbl => Some(&self.xor_mbl),
						Step::Tybox if self.uses_xor => Some(&self.xor),
						_ => None,
					};
					if let Some(xor) = xor {
						let xor = xor.partial_map(r, col);

						let n01 = |v: SColumn, n: HighLow| {
							let a = xor.map(
								Purpose::High,
								v,
								n,
								aa.column_nibble(v, n),
								bb.column_nibble(v, n),
							);
							let b = xor.map(
								Purpose::Low,
								v,
								n,
								cc.column_nibble(v, n),
								dd.column_nibble(v, n),
							);

							xor.map(Purpose::Output, v, n, a, b)
						};

						let n0123 =
							|v: SColumn| X::nibs(n01(v, HighLow::High), n01(v, HighLow::Low));

						data.set_column(col, Column(SColumn::ALL.map(n0123)));
					} else {
						let n0123 = |v: SColumn| aa[v] ^ bb[v] ^ cc[v] ^ dd[v];
						data.set_column(col, Column(SColumn::ALL.map(n0123)));
					}
				}
			}
		}
		shift_rows(
			data,
			if self.inv {
				&INV_SHIFT_ROWS_TAB
			} else {
				&SHIFT_ROWS_TAB
			},
		);
		// sub_bytes + add_round_key
		self.tboxes_last.apply(data);
	}
}

#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct State(StateMap<X>);
impl State {
	pub const TWO_ONE_NINE_TWO_TEST_VECTOR: Self = State::from_bytes(*b"Two One Nine Two");
	pub const TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR: Self =
		State::from_bytes(*b"\x29\xc3\x50\x5f\x57\x14\x20\xf6\x40\x22\x99\xb3\x1a\x02\xd7\x3a");
	pub const TWO_ONE_NINE_TWO_AES256_KUNG_FU_TEST_VECTOR: Self =
		State::from_bytes(*b"\xe3\x5a\x6d\xcb\x19\xb2\x01\xa0\x1e\xbc\xfa\x8a\xa2\x2b\x57\x59");

	const fn from_bytes(data: [u8; 16]) -> Self {
		Self(unsafe { transmute::<[u8; 16], StateMap<X>>(data) })
	}

	fn from_x(x: X) -> Self {
		Self(StateMap([x; 16]))
	}

	fn set_column(&mut self, r: SColumn, v: Column) {
		for row in SRow::all() {
			self[SPos::column_row(r, row)] = v[row];
		}
	}
	fn get_column(&self, r: SColumn) -> Column {
		Column(SRow::ALL.map(|c| self[SPos::column_row(r, c)]))
	}
}
impl Index<SPos> for State {
	type Output = X;

	fn index(&self, index: SPos) -> &Self::Output {
		&self.0[index]
	}
}
impl IndexMut<SPos> for State {
	fn index_mut(&mut self, index: SPos) -> &mut Self::Output {
		&mut self.0[index]
	}
}

macro_rules! u2_newtype {
	($id:ident) => {
		#[derive(Clone, Copy, Default)]
		pub struct $id(U2);

		impl $id {
			const ALL: [Self; 4] = [Self::_0, Self::_1, Self::_2, Self::_3];
			const _0: Self = Self(U2::_0);
			const _1: Self = Self(U2::_1);
			const _2: Self = Self(U2::_2);
			const _3: Self = Self(U2::_3);
			fn all() -> impl Iterator<Item = Self> {
				U2::all().map($id)
			}
			fn as_index(&self) -> usize {
				self.0.as_index()
			}
			fn from_index(i: usize) -> Self {
				Self::ALL[i]
			}
		}
	};
}

u2_newtype!(SColumn);
u2_newtype!(SRow);

#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub struct SPos(U4);
impl SPos {
	const ALL: [Self; 16] = [
		Self(U4::_0),
		Self(U4::_1),
		Self(U4::_2),
		Self(U4::_3),
		Self(U4::_4),
		Self(U4::_5),
		Self(U4::_6),
		Self(U4::_7),
		Self(U4::_8),
		Self(U4::_9),
		Self(U4::_10),
		Self(U4::_11),
		Self(U4::_12),
		Self(U4::_13),
		Self(U4::_14),
		Self(U4::_15),
	];
	fn column_row(column: SColumn, row: SRow) -> Self {
		Self(U4::ji(column.0, row.0))
	}
	fn all() -> impl Iterator<Item = Self> {
		U4::all().map(Self)
	}
	fn row(&self) -> SRow {
		SRow(self.0.i())
	}
	fn column(&self) -> SColumn {
		SColumn(self.0.j())
	}
	fn as_index(&self) -> usize {
		self.0.as_index()
	}
	fn from_index(v: usize) -> Self {
		Self(U4::from_index(v))
	}
}

// 0 - 4
#[derive(Clone, Copy, Default)]
#[repr(u8)]
#[derive(Debug, PartialEq)]
pub enum U2 {
	#[default]
	_0 = 0,
	_1 = 1,
	_2 = 2,
	_3 = 3,
}
impl U2 {
	const ALL: [Self; 4] = [Self::_0, Self::_1, Self::_2, Self::_3];
	fn all() -> impl Iterator<Item = Self> {
		Self::ALL.into_iter()
	}

	const fn as_index(self) -> usize {
		self as usize
	}
	pub const fn from_index(index: usize) -> Self {
		Self::ALL[index]
	}
}
impl BitXor for U2 {
	type Output = U2;

	fn bitxor(self, rhs: Self) -> Self::Output {
		Self::from_index(self.as_index() ^ rhs.as_index())
	}
}

#[derive(Clone, Copy, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum U4 {
	#[default]
	_0 = 0,
	_1 = 1,
	_2 = 2,
	_3 = 3,
	_4 = 4,
	_5 = 5,
	_6 = 6,
	_7 = 7,
	_8 = 8,
	_9 = 9,
	_10 = 10,
	_11 = 11,
	_12 = 12,
	_13 = 13,
	_14 = 14,
	_15 = 15,
}
impl U4 {
	const ALL: [Self; 16] = [
		Self::_0,
		Self::_1,
		Self::_2,
		Self::_3,
		Self::_4,
		Self::_5,
		Self::_6,
		Self::_7,
		Self::_8,
		Self::_9,
		Self::_10,
		Self::_11,
		Self::_12,
		Self::_13,
		Self::_14,
		Self::_15,
	];
	const fn ji(j: U2, i: U2) -> Self {
		Self::ALL[(j.as_index() << 2) | i.as_index()]
	}
	fn j(self) -> U2 {
		let i = self.as_index();
		U2::ALL[i >> 2]
	}
	fn i(self) -> U2 {
		let i = self.as_index();
		U2::ALL[i & 0b11]
	}

	fn all() -> impl Iterator<Item = Self> {
		Self::ALL.iter().copied()
	}

	pub const fn as_index(self) -> usize {
		self as usize
	}
	pub const fn from_index(v: usize) -> Self {
		Self::ALL[v]
	}
}
impl fmt::Debug for U4 {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{:x}", self.as_index())
	}
}
impl BitXor for U4 {
	type Output = U4;

	fn bitxor(self, rhs: Self) -> Self::Output {
		Self::from_index(self as usize ^ rhs as usize)
	}
}

// Performs the AddRoundKey step. Each round has its own pre-generated 16-byte key in the
// form of 4 integers (the "w" array). Each integer is XOR'd by one column of the state.
// Also performs the job of InvAddRoundKey(); since the function is a simple XOR process,
// it is its own inverse.
pub fn add_round_key(state: &mut State, round_key: &RoundKey) {
	add_shifted_round_key(state, round_key, &ShiftRowsBijection::IDENTITY)
}
fn add_shifted_round_key(state: &mut State, round_key: &RoundKey, shift: &ShiftRowsBijection) {
	for pos in SPos::all() {
		state[pos] ^= round_key.0[shift.map(pos).0];
	}
}

fn sub_bytes<S: SBox>(sbox: S, state: &mut State) {
	for (_, ele) in state.0.iter_mut() {
		*ele = sbox.sub_byte(*ele)
	}
}
fn gf_mul(a: X, b: u8, poly: dual::Poly) -> X {
	X(consts::gf_mul_pol_slow(a.0, b, poly))
}

pub fn mix_columns(state: &mut State, coeffs: &dual::MixColCoeffs, poly: dual::Poly) {
	for column in SColumn::all() {
		let Column([a, b, c, d]) = state.get_column(column);

		state.set_column(
			column,
			Column([
				gf_mul(a, coeffs.c2, poly) ^ gf_mul(b, coeffs.c3, poly) ^ c ^ d,
				a ^ gf_mul(b, coeffs.c2, poly) ^ gf_mul(c, coeffs.c3, poly) ^ d,
				a ^ b ^ gf_mul(c, coeffs.c2, poly) ^ gf_mul(d, coeffs.c3, poly),
				gf_mul(a, coeffs.c3, poly) ^ b ^ c ^ gf_mul(d, coeffs.c2, poly),
			]),
		);
	}
}

pub fn mix_columns_inv(state: &mut State, coeffs: &dual::MixColCoeffs, poly: dual::Poly) {
	for column in SColumn::all() {
		let Column([a, b, c, d]) = state.get_column(column);

		state.set_column(
			column,
			Column([
				gf_mul(a, coeffs.c14, poly)
					^ gf_mul(b, coeffs.c11, poly)
					^ gf_mul(c, coeffs.c13, poly)
					^ gf_mul(d, coeffs.c9, poly),
				gf_mul(a, coeffs.c9, poly)
					^ gf_mul(b, coeffs.c14, poly)
					^ gf_mul(c, coeffs.c11, poly)
					^ gf_mul(d, coeffs.c13, poly),
				gf_mul(a, coeffs.c13, poly)
					^ gf_mul(b, coeffs.c9, poly)
					^ gf_mul(c, coeffs.c14, poly)
					^ gf_mul(d, coeffs.c11, poly),
				gf_mul(a, coeffs.c11, poly)
					^ gf_mul(b, coeffs.c13, poly)
					^ gf_mul(c, coeffs.c9, poly)
					^ gf_mul(d, coeffs.c14, poly),
			]),
		);
	}
}

/// Round index
#[derive(PartialEq, Eq, Clone, Copy, PartialOrd, Ord, Debug)]
pub struct RI(usize);
impl RI {
	fn all<const RN: usize>() -> impl Iterator<Item = RI> {
		(0..RN).map(RI)
	}
}

#[derive(Clone, Copy)]
#[repr(u8)]
pub enum Purpose {
	High = 0,
	Low = 1,
	Output = 2,
}
impl Purpose {
	const ALL: [Self; 3] = [Purpose::High, Purpose::Low, Purpose::Output];
	fn all() -> impl Iterator<Item = Self> {
		Self::ALL.into_iter()
	}
	fn as_index(self) -> usize {
		self as usize
	}
	fn from_index(i: usize) -> Self {
		Self::ALL[i]
	}
}

#[derive(Clone, Copy)]
#[repr(u8)]
pub enum HighLow {
	High = 0,
	Low = 1,
}
impl HighLow {
	const ALL: [Self; 2] = [Self::Low, Self::High];
	fn all() -> impl Iterator<Item = Self> {
		Self::ALL.into_iter()
	}
	fn as_index(self) -> usize {
		self as usize
	}
	fn from_index(i: usize) -> Self {
		Self::ALL[i]
	}
}
impl<T, const N: usize> Distribution<RIArr<T, N>> for StandardUniform
where
	StandardUniform: Distribution<T>,
{
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> RIArr<T, N> {
		let mut out = [const { MaybeUninit::<T>::uninit() }; N];
		for ele in out.iter_mut() {
			ele.write(rng.random());
		}
		RIArr(unsafe { std::mem::transmute_copy(&out) })
	}
}

pub fn try_array_from_fn<const S: usize, I, T, E>(
	from_index: impl Fn(usize) -> I,
	mut f: impl FnMut(I) -> Result<T, E>,
) -> Result<[T; S], E> {
	let mut out = [const { <MaybeUninit<T>>::uninit() }; S];
	let mut initialized = 0;
	for (i, v) in out.iter_mut().enumerate() {
		let new_v = f(from_index(i));
		match new_v {
			Ok(new_v) => {
				v.write(new_v);
				initialized += 1;
			}
			Err(e) => {
				for i in 0..initialized {
					unsafe { out[i].assume_init_drop() }
				}
				return Err(e);
			}
		}
	}
	Ok(unsafe { std::mem::transmute_copy::<_, [T; S]>(&out) })
}
macro_rules! fixed_map {
	($id:ident($s:literal as $t:ty)) => {
		#[derive(Clone, Copy, Debug, PartialEq)]
		pub struct $id<T>([T; $s]);
		impl<T> $id<T> {
			pub fn from_fn(f: impl FnMut($t) -> T) -> Self {
				Self(<$t>::ALL.map(f))
			}
			pub fn try_from_fn<E>(f: impl FnMut($t) -> Result<T, E>) -> Result<Self, E> {
				try_array_from_fn(<$t>::from_index, f).map(Self)
			}
			pub fn iter_mut(&mut self) -> impl Iterator<Item = ($t, &mut T)> {
				<$t>::all().zip(self.0.iter_mut())
			}
		}
		impl<T> Index<$t> for $id<T> {
			type Output = T;
			fn index(&self, index: $t) -> &Self::Output {
				&self.0[index.as_index()]
			}
		}
		impl<T> IndexMut<$t> for $id<T> {
			fn index_mut(&mut self, index: $t) -> &mut Self::Output {
				&mut self.0[index.as_index()]
			}
		}
		impl<T> Default for $id<T>
		where
			T: Default,
		{
			fn default() -> Self {
				Self(array::from_fn(|_| Default::default()))
			}
		}
		impl<T> Distribution<$id<T>> for StandardUniform
		where
			StandardUniform: Distribution<T>,
			T: Default,
		{
			fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> $id<T> {
				$id(std::array::from_fn(|_| rng.random()))
			}
		}
	};
	($id:ident<$v:ty>($s:literal as $t:ty)) => {
		#[derive(Clone, Copy, PartialEq)]
		pub struct $id([$v; $s]);
		impl $id {
			pub fn iter_mut(&mut self) -> impl Iterator<Item = ($t, &mut $v)> {
				<$t>::all().zip(self.0.iter_mut())
			}
		}
		impl Index<$t> for $id {
			type Output = $v;
			fn index(&self, index: $t) -> &Self::Output {
				&self.0[index.as_index()]
			}
		}
		impl IndexMut<$t> for $id {
			fn index_mut(&mut self, index: $t) -> &mut Self::Output {
				&mut self.0[index.as_index()]
			}
		}
		impl Default for $id {
			fn default() -> Self {
				Self([Default::default(); $s])
			}
		}
	};
}
fixed_map!(NibbleMap(16 as U4));
fixed_map!(StateMap(16 as SPos));
fixed_map!(DibitMap(4 as U2));
fixed_map!(RowMap(4 as SRow));
fixed_map!(ColumnMap(4 as SColumn));
fixed_map!(XArr(256 as X));
fixed_map!(PurposeMap(3 as Purpose));
fixed_map!(HighLowMap(2 as HighLow));

/// Array of 256 elements
impl XArr<X> {
	pub const fn new_static(v: [u8; 256]) -> Self {
		Self(unsafe { std::mem::transmute::<[u8; 256], [X; 256]>(v) })
	}
	pub fn from_fn_par<F>(map: F) -> Self
	where
		F: Fn(X) -> X,
		F: Sync,
	{
		let mut out = [X(0); 256];
		out.par_iter_mut()
			.enumerate()
			.for_each(|(i, p)| *p = map(X(i as u8)));
		Self(out)
	}
}

/// State byte, GF(2^256) value
#[derive(Clone, Copy, Default, PartialEq)]
pub struct X(u8);
impl X {
	const ALL: [Self; 256] = {
		let mut out = [X(0); 256];
		let mut i = 0;
		while i <= 255 {
			out[i] = X(i as u8);
			i += 1;
		}
		out
	};
	fn as_index(&self) -> usize {
		self.0 as usize
	}
	fn from_index(i: usize) -> Self {
		Self::ALL[i]
	}
	fn nibs(a: U4, b: U4) -> Self {
		Self(((a.as_index() << 4) | b.as_index()) as u8)
	}
	fn as_nibs(self) -> (U4, U4) {
		let a = self.0 >> 4;
		let b = self.0 & 0xf;
		(U4::from_index(a as usize), U4::from_index(b as usize))
	}
	fn all() -> impl Iterator<Item = X> {
		(0..=255).map(X)
	}
}
impl BitXor for X {
	type Output = X;

	fn bitxor(self, rhs: Self) -> Self::Output {
		X(self.0 ^ rhs.0)
	}
}
impl BitXorAssign for X {
	fn bitxor_assign(&mut self, rhs: Self) {
		*self = *self ^ rhs;
	}
}
impl fmt::Debug for X {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{:0>2x}", self.0)
	}
}

/// Array of [number of rounds - 1] elements
#[derive(Debug)]
pub struct RIArr<T, const NRM1: usize>(pub(crate) [T; NRM1]);
impl<T: Clone, const NRM1: usize> Clone for RIArr<T, NRM1> {
	fn clone(&self) -> Self {
		Self(self.0.clone())
	}
}
impl<T, const NRM1: usize> RIArr<T, NRM1> {
	fn from_fn(mut f: impl FnMut(RI) -> T) -> Self {
		Self(array::from_fn(|ri| f(RI(ri))))
	}

	fn iter_mut(&mut self) -> impl DoubleEndedIterator<Item = (RI, &mut T)> + ExactSizeIterator {
		self.0.iter_mut().enumerate().map(|(ri, v)| (RI(ri), v))
	}
	pub fn try_from_fn<E>(f: impl FnMut(RI) -> Result<T, E>) -> Result<Self, E> {
		try_array_from_fn(RI, f).map(Self)
	}
}
impl<T, const NRM1: usize> Index<RI> for RIArr<T, NRM1> {
	type Output = T;

	fn index(&self, index: RI) -> &Self::Output {
		&self.0[index.0]
	}
}
impl<T, const NRM1: usize> IndexMut<RI> for RIArr<T, NRM1> {
	fn index_mut(&mut self, index: RI) -> &mut Self::Output {
		&mut self.0[index.0]
	}
}
impl<T, const NRM1: usize> Default for RIArr<T, NRM1>
where
	T: Default + Copy,
{
	fn default() -> Self {
		Self([T::default(); NRM1])
	}
}

fixed_map!(Row<X>(4 as SColumn));
fixed_map!(Column<X>(4 as SRow));

impl fmt::Debug for Row {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "Row(")?;
		for v in self.0 {
			write!(f, "{:0>2x}", v.0)?;
		}
		write!(f, ")")
	}
}

impl Row {
	fn column_nibble(&self, column: SColumn, high_low: HighLow) -> U4 {
		let v = self[column];
		let (high, low) = v.as_nibs();
		match high_low {
			HighLow::High => high,
			HighLow::Low => low,
		}
	}
	fn set_column_nibble(&mut self, column: SColumn, high_low: HighLow, to: U4) {
		let v = self[column];
		let (mut high, mut low) = v.as_nibs();
		match high_low {
			HighLow::High => high = to,
			HighLow::Low => low = to,
		};
		self[column] = X::nibs(high, low)
	}
	fn as_bytes(&self) -> [u8; 4] {
		[self.0[0].0, self.0[1].0, self.0[2].0, self.0[3].0]
	}
	fn from_bytes(v: [u8; 4]) -> Self {
		Self([X(v[0]), X(v[1]), X(v[2]), X(v[3])])
	}
}
impl Column {
	fn as_bytes(&self) -> [u8; 4] {
		[self.0[0].0, self.0[1].0, self.0[2].0, self.0[3].0]
	}
	fn from_bytes(v: [u8; 4]) -> Self {
		Self([X(v[0]), X(v[1]), X(v[2]), X(v[3])])
	}
}

fn rot_word(word: Word) -> Word {
	let [a, b, c, d] = word.to_bytes();
	Word::from_bytes([b, c, d, a])
}
fn sub_word<S: SBox>(sbox: S, word: Word) -> Word {
	Word::from_bytes(word.to_bytes().map(|b| sbox.sub_byte(X(b)).0))
}

#[derive(Copy, Clone, Zeroize)]
#[repr(transparent)]
pub struct Word([u8; 4]);
impl Word {
	fn from_bytes(b: [u8; 4]) -> Self {
		Self(b)
	}
	fn to_bytes(self) -> [u8; 4] {
		self.0
	}
}
impl BitXor for Word {
	type Output = Word;

	fn bitxor(self, rhs: Self) -> Self::Output {
		let mut b = self.to_bytes();
		let c = rhs.to_bytes();
		for i in 0..4 {
			b[i] ^= c[i];
		}
		Self::from_bytes(b)
	}
}

macro_rules! impl_aes_tables {
	($name:ident, $keys:ident, $nrm1:literal, $key:ident) => {
		pub type $name = Tables<$nrm1>;
		pub type $keys = RoundKeys<$nrm1>;
		impl $name {
			pub const ROUNDS: usize = $nrm1;
			pub fn from_key(key: $key, inv: bool) -> Self {
				Self::from_nonstandard_key(
					&PrecomputedSBox::aes_standard(),
					key,
					inv,
					Dual::STANDARD,
				)
			}
			pub fn from_unchecked_round_keys(round_keys: &$keys, inv: bool) -> Self {
				Self::from_nonstandard_round_keys(
					&PrecomputedSBox::aes_standard(),
					round_keys,
					inv,
					Dual::STANDARD,
				)
			}
		}
	};
}

impl_aes_tables!(Aes128Tables, Aes128RoundKeys, 9, Aes128Key);
impl_aes_tables!(Aes192Tables, Aes192RoundKeys, 11, Aes192Key);
impl_aes_tables!(Aes256Tables, Aes256RoundKeys, 13, Aes256Key);

#[derive(Clone, Copy)]
pub struct RoundKey(NibbleMap<X>);
impl RoundKey {
	fn get_word(&self, j: U2) -> Word {
		Word(U2::ALL.map(|i| self.0[U4::ji(j, i)].0))
	}
	fn set_word(&mut self, i: U2, v: Word) {
		let [a, b, c, d] = v.to_bytes();

		self.0[U4::ji(i, U2::_0)] = X(a);
		self.0[U4::ji(i, U2::_1)] = X(b);
		self.0[U4::ji(i, U2::_2)] = X(c);
		self.0[U4::ji(i, U2::_3)] = X(d)
	}
	fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for ele in self.0.0 {
			out.write_all(&[ele.0])?
		}
		Ok(())
	}
	fn read_from(mut input: impl Read) -> io::Result<Self> {
		Ok(Self(NibbleMap::try_from_fn(|_| {
			let mut wdata = [0];
			input.read_exact(&mut wdata)?;
			Ok(X(wdata[0])) as io::Result<_>
		})?))
	}
}
impl Distribution<RoundKey> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> RoundKey {
		RoundKey(NibbleMap::from_fn(|_| X(rng.random())))
	}
}
impl fmt::Debug for RoundKey {
	fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(fmt, "RoundKey(")?;
		for v in self.0.0 {
			write!(fmt, "{v:?}")?;
		}
		write!(fmt, ")")
	}
}

pub fn encrypt<const RNM1: usize>(state: &mut State, rk: &RoundKeys<RNM1>) {
	encrypt_nonstandard(&PrecomputedSBox::aes_standard(), state, rk, Dual::STANDARD)
}
pub fn encrypt_nonstandard<S: SBox + Copy, const RNM1: usize>(
	sbox: S,
	state: &mut State,
	rk: &RoundKeys<RNM1>,
	dual: Dual,
) {
	let coeffs = MixColCoeffs::for_dual(dual);
	let poly = dual.poly;

	let q = Q::for_dual(dual);
	state.apply_q(q);

	add_round_key(state, &rk.round(RI(0)));

	for r in 1..=RNM1 {
		let r = RI(r);
		sub_bytes(sbox, state);
		shift_rows(state, &SHIFT_ROWS_TAB);
		mix_columns(state, &coeffs, poly);
		add_round_key(state, &rk.round(r))
	}

	sub_bytes(sbox, state);
	shift_rows(state, &SHIFT_ROWS_TAB);
	add_round_key(state, &rk.round(RI(RNM1 + 1)));
}
pub fn decrypt<const RNM1: usize>(state: &mut State, rk: &RoundKeys<RNM1>) {
	decrypt_nonstandard(&PrecomputedSBox::aes_standard(), state, rk, Dual::STANDARD);
}
pub fn decrypt_nonstandard<S: SBox + Copy, const RNM1: usize>(
	sbox: S,
	state: &mut State,
	rk: &RoundKeys<RNM1>,
	config: Dual,
) {
	let inv_sbox = PrecomputedSBox::inversed(&sbox);
	let coeffs = MixColCoeffs::for_dual(config);
	let poly = config.poly;

	add_round_key(state, &rk.round(RI(RNM1 + 1)));

	for r in (1..=RNM1).rev() {
		let r = RI(r);
		shift_rows(state, &INV_SHIFT_ROWS_TAB);
		sub_bytes(&inv_sbox, state);
		add_round_key(state, &rk.round(r));
		mix_columns_inv(state, &coeffs, poly);
	}

	shift_rows(state, &INV_SHIFT_ROWS_TAB);
	sub_bytes(inv_sbox, state);
	add_round_key(state, &rk.round(RI(0)));

	let q = Q::for_dual(config);
	state.apply_q_inv(q);
}

#[cfg(test)]
mod tests {
	use rand::Rng;
	use rand::rng;
	use tracing::debug;

	use crate::Aes256Tables;
	use crate::Security;
	use crate::decrypt;
	use crate::encoding::ExternalEncoding;
	use crate::encoding::apply_encoding;
	use crate::encrypt;
	use crate::hardware::decrypt_hardware;
	use crate::hardware::encrypt_hardware;
	use crate::key::Aes128Key;
	use crate::key::Aes256Key;
	use crate::{Aes128Tables, IV, State};

	#[test_log::test]
	fn normal() {
		let mut tables = Aes128Tables::from_key(Aes128Key::NIST_CFB_E2, false);
		tables.apply_security(Security::full(), &mut rng());

		const NIST_MESSAGE_CFB_E2: [u8; 64] = [
			0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93,
			0x17, 0x2a, 0xae, 0x2d, 0x8a, 0x57, 0x1e, 0x03, 0xac, 0x9c, 0x9e, 0xb7, 0x6f, 0xac,
			0x45, 0xaf, 0x8e, 0x51, 0x30, 0xc8, 0x1c, 0x46, 0xa3, 0x5c, 0xe4, 0x11, 0xe5, 0xfb,
			0xc1, 0x19, 0x1a, 0x0a, 0x52, 0xef, 0xf6, 0x9f, 0x24, 0x45, 0xdf, 0x4f, 0x9b, 0x17,
			0xad, 0x2b, 0x41, 0x7b, 0xe6, 0x6c, 0x37, 0x10,
		];

		const NIST_TEST_RESULT_CFB_E2: [u8; 64] = [
			0x3b, 0x3f, 0xd9, 0x2e, 0xb7, 0x2d, 0xad, 0x20, 0x33, 0x34, 0x49, 0xf8, 0xe8, 0x3c,
			0xfb, 0x4a, 0xc8, 0xa6, 0x45, 0x37, 0xa0, 0xb3, 0xa9, 0x3f, 0xcd, 0xe3, 0xcd, 0xad,
			0x9f, 0x1c, 0xe5, 0x8b, 0x26, 0x75, 0x1f, 0x67, 0xa3, 0xcb, 0xb1, 0x40, 0xb1, 0x80,
			0x8c, 0xf1, 0x87, 0xa4, 0xf4, 0xdf, 0xc0, 0x4b, 0x05, 0x35, 0x7c, 0x5d, 0x1c, 0x0e,
			0xea, 0xc4, 0xc6, 0x6f, 0x9f, 0xf7, 0xf2, 0xe6,
		];

		let mut message = NIST_MESSAGE_CFB_E2;
		tables.encrypt_cfb(IV::NIST_CFB_E1, &mut message);
		debug!("encrypted = {message:?}");
		assert_eq!(
			message, NIST_TEST_RESULT_CFB_E2,
			"encryption is unsuccessful"
		);
		tables.decrypt_cfb(IV::NIST_CFB_E1, &mut message);
		assert_eq!(message, NIST_MESSAGE_CFB_E2, "decryption is unsuccessful");
	}

	#[test]
	fn inv() {
		let s = Security::full();
		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, false);
		tables.apply_security(s, &mut rng());
		let mut tables_inv = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, true);
		tables_inv.apply_security(s, &mut rng());

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
		tables_inv.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR)
	}

	#[test]
	fn external_encoding() {
		let rng = &mut rng();
		let s = Security::full();

		let encoding: ExternalEncoding = rng.random();
		let encoding2: ExternalEncoding = rng.random();

		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, false);
		tables.apply_security(s, rng);
		tables.output_encoding(&encoding, false);
		tables.output_encoding(&encoding2, false);

		let mut tables_inv = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, true);
		tables_inv.apply_security(s, rng);
		tables_inv.input_encoding(&encoding, true);
		tables_inv.input_encoding(&encoding2, true);

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_ne!(
			data,
			State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR,
			"encoded aes output should not match original aes"
		);
		tables_inv.cipher(&mut data);
		assert_eq!(
			data,
			State::TWO_ONE_NINE_TWO_TEST_VECTOR,
			"with encodings the values should roundtrip"
		);
	}

	#[test]
	fn external_encoding_standalone() {
		let rng = &mut rng();
		let s = Security::full();

		let encoding: ExternalEncoding = rng.random();

		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, false);
		tables.apply_security(s, rng);
		tables.input_encoding(&encoding, true);

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		apply_encoding(&mut data, &encoding, false);

		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR,);
	}
	#[test]
	fn external_decoding_standalone() {
		let rng = &mut rng();
		let s = Security::full();

		let encoding: ExternalEncoding = rng.random();

		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, true);
		tables.apply_security(s, rng);
		tables.output_encoding(&encoding, false);

		let mut data = State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR;
		tables.cipher(&mut data);
		apply_encoding(&mut data, &encoding, true);

		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR,);
	}

	#[test]
	fn self_cancelling_external_decoding_standalone() {
		let rng = &mut rng();
		let s = Security::full();

		let input_encoding: ExternalEncoding = rng.random();
		let output_encoding: ExternalEncoding = rng.random();

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		// decode(encode(v)) is isomorphic to encode(decode(v))
		apply_encoding(&mut data, &output_encoding, true);

		let round_keys = Aes128Key::KUNG_FU_TEST_VECTOR.expand();
		encrypt(&mut data, &round_keys);

		apply_encoding(&mut data, &input_encoding, false);

		let mut tables = Aes128Tables::from_unchecked_round_keys(&round_keys, true);
		tables.apply_security(s, rng);
		tables.input_encoding(&input_encoding, true);
		tables.output_encoding(&output_encoding, false);

		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR,);
	}

	#[test]
	fn inv_256() {
		let mut key = [0; 32];
		key[0] = 0x80;
		let s = Security::full();
		let mut tables = Aes256Tables::from_key(Aes256Key::new(key), false);
		tables.apply_security(s, &mut rng());
		let mut tables_inv = Aes256Tables::from_key(Aes256Key::new(key), true);
		tables_inv.apply_security(s, &mut rng());

		let mut data = State::from_bytes([0; 16]);
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES256_KUNG_FU_TEST_VECTOR);
		tables_inv.cipher(&mut data);
		assert_eq!(data, State::from_bytes([0; 16]));
	}

	#[test]
	fn internal_encodings() {
		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, false);
		tables.apply_security(Security::full(), &mut rng());

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
	}
	#[test]
	fn internal_encodings_inv() {
		let mut tables = Aes128Tables::from_key(Aes128Key::KUNG_FU_TEST_VECTOR, true);
		tables.apply_security(Security::full(), &mut rng());

		let mut data = State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR;
		tables.cipher(&mut data);
		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR);
	}

	#[test_log::test]
	fn standard() {
		let round_keys = Aes128Key::KUNG_FU_TEST_VECTOR.expand();

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		encrypt(&mut data, &round_keys);

		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR,);

		decrypt(&mut data, &round_keys);

		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR);
	}
	#[test]
	fn hardware() {
		let round_keys = Aes128Key::KUNG_FU_TEST_VECTOR.expand();

		let mut data = State::TWO_ONE_NINE_TWO_TEST_VECTOR;
		assert!(
			std::arch::is_x86_feature_detected!("sse2")
				&& std::arch::is_x86_feature_detected!("aes")
		);
		// SAFETY: Asserted that those intrinsics are supported
		unsafe {
			encrypt_hardware(&mut data, &round_keys);
		};

		assert_eq!(data, State::TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR,);

		unsafe { decrypt_hardware(&mut data, &round_keys) }

		assert_eq!(data, State::TWO_ONE_NINE_TWO_TEST_VECTOR,);
	}
}
