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

use crate::key::RoundKeys;
use crate::{RI, State, add_round_key};

/// # SAFETY
///
/// This functions uses AES and SSE2 intrinsics, they should be available
#[target_feature(enable = "aes")]
#[target_feature(enable = "sse2")]
pub unsafe fn encrypt_hardware<const RNM1: usize>(state: &mut State, rk: &RoundKeys<RNM1>) {
	use std::arch::x86_64::{
		_mm_aesenc_si128, _mm_aesenclast_si128, _mm_loadu_si128, _mm_storeu_si128,
	};

	add_round_key(state, &rk.round(RI(0)));

	let mut s = unsafe { _mm_loadu_si128(state.0.0.as_ptr().cast()) };

	for r in 1..=RNM1 {
		let r = RI(r);

		let k = unsafe { _mm_loadu_si128(rk.round(r).0.0.as_ptr().cast()) };
		s = _mm_aesenc_si128(s, k);
	}

	let kf = unsafe { _mm_loadu_si128(rk.round(RI(RNM1 + 1)).0.0.as_ptr().cast()) };
	s = _mm_aesenclast_si128(s, kf);

	unsafe { _mm_storeu_si128(state.0.0.as_mut_ptr().cast(), s) };
}
/// # SAFETY
///
/// This functions uses AES and SSE2 intrinsics, they should be available
#[target_feature(enable = "aes")]
#[target_feature(enable = "sse2")]
pub unsafe fn decrypt_hardware<const RNM1: usize>(state: &mut State, rk: &RoundKeys<RNM1>) {
	use std::arch::x86_64::{
		_mm_aesdec_si128, _mm_aesdeclast_si128, _mm_aesimc_si128, _mm_loadu_si128, _mm_storeu_si128,
	};

	add_round_key(state, &rk.round(RI(RNM1 + 1)));

	let mut s = unsafe { _mm_loadu_si128(state.0.0.as_ptr().cast()) };

	for r in (1..=RNM1).rev() {
		let r = RI(r);

		// shift_rows(inv)
		// sub_bytes(inv_sbox)
		// add_round_key
		// mix_columns_inv
		let k = unsafe { _mm_loadu_si128(rk.round(r).0.0.as_ptr().cast()) };
		// Inverse mix columns needs to be applied on key before decryption.
		// TODO: Use functions from this crate instead of hardware-accelerated one.
		let k = _mm_aesimc_si128(k);
		s = _mm_aesdec_si128(s, k);
	}

	// shift_rows(inv)
	// sub_bytes(inv_sbox)
	// add_round_key
	let kf = unsafe { _mm_loadu_si128(rk.round(RI(0)).0.0.as_ptr().cast()) };
	s = _mm_aesdeclast_si128(s, kf);

	unsafe { _mm_storeu_si128(state.0.0.as_mut_ptr().cast(), s) };
}
