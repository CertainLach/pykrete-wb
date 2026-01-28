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
use std::mem::transmute;

use num_integer::Integer as _;
use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use tracing::{Level, debug, instrument};
use zeroize::ZeroizeOnDrop;

use crate::sbox::{PrecomputedSBox, SBox};
use crate::{NibbleMap, RI, RIArr, RoundKey, U2, Word, X, rot_word, sub_word};

#[derive(ZeroizeOnDrop)]
pub struct Key<const NK: usize>(pub(crate) [Word; NK]);

macro_rules! impl_aes_key {
	($name:ident, $words:literal, $bytes:literal, $nrm1:literal) => {
		pub type $name = Key<$words>;
		impl $name {
			pub const NRM1: usize = $nrm1;

			pub const fn new(v: [u8; $bytes]) -> Self {
				Self(unsafe { transmute::<[u8; $bytes], [Word; $words]>(v) })
			}

			pub fn expand(self) -> RoundKeys<$nrm1> {
				expand_nonstandard_keys(
					&PrecomputedSBox::aes_standard(),
					self,
					Dual::STANDARD,
				)
			}
		}
	};
}

impl_aes_key!(Aes128Key, 4, 16, 9);
impl_aes_key!(Aes192Key, 6, 24, 11);
impl_aes_key!(Aes256Key, 8, 32, 13);

impl Aes128Key {
	#[cfg(test)]
	pub const NIST_CFB_E1: Self = Self::new([
		0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f,
		0x3c,
	]);
	pub const NIST_CFB_E2: Self = Self::new([
		0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f,
		0x3c,
	]);
	pub const KUNG_FU_TEST_VECTOR: Self = Self::new(*b"Thats my Kung Fu");
}

impl<const NK: usize> Key<NK> {
	pub fn word(&self, round: usize) -> Word {
		self.0[round]
	}
	pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for ele in self.0 {
			out.write_all(&ele.0)?;
		}
		Ok(())
	}
	pub fn read_from(mut input: impl Read) -> io::Result<Self> {
		let mut out = [Word([0, 0, 0, 0]); NK];
		for ele in out.iter_mut() {
			let mut wdata = [0, 0, 0, 0];
			input.read_exact(&mut wdata)?;
			*ele = Word(wdata)
		}
		Ok(Self(out))
	}
}

#[derive(Debug)]
pub struct RoundKeys<const NRM1: usize> {
	pub rounds: RIArr<RoundKey, NRM1>,
	pub last: (RoundKey, RoundKey),
}
impl<const NRM1: usize> RoundKeys<NRM1> {
	pub(crate) fn new() -> Self {
		Self {
			rounds: RIArr([const { RoundKey(NibbleMap([X(0); 16])) }; NRM1]),
			last: (
				RoundKey(NibbleMap([X(0); 16])),
				RoundKey(NibbleMap([X(0); 16])),
			),
		}
	}
	pub fn round_mut(&mut self, i: RI) -> &mut RoundKey {
		if i == RI(NRM1) {
			&mut self.last.0
		} else if i == RI(NRM1 + 1) {
			&mut self.last.1
		} else {
			&mut self.rounds[i]
		}
	}
	pub fn round(&self, i: RI) -> RoundKey {
		if i.0 == NRM1 {
			self.last.0
		} else if i.0 == NRM1 + 1 {
			self.last.1
		} else {
			self.rounds[i]
		}
	}
	pub fn get_schedule_word(&self, i: usize) -> Word {
		let (round, word) = i.div_mod_floor(&4);
		self.round(RI(round)).get_word(U2::from_index(word))
	}
	pub fn set_schedule_word(&mut self, i: usize, v: Word) {
		let (round, word) = i.div_mod_floor(&4);
		self.round_mut(RI(round)).set_word(U2::from_index(word), v)
	}

	pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
		for round in &self.rounds.0 {
			round.write_to(&mut out)?;
		}
		let (a, b) = &self.last;
		a.write_to(&mut out)?;
		b.write_to(out)
	}
	pub fn read_from(mut input: impl Read) -> io::Result<Self> {
		let rounds = RIArr::try_from_fn(|_| RoundKey::read_from(&mut input))?;
		let a = RoundKey::read_from(&mut input)?;
		let b = RoundKey::read_from(input)?;
		Ok(Self {
			rounds,
			last: (a, b),
		})
	}
}
impl<const NRM1: usize> Distribution<RoundKeys<NRM1>> for StandardUniform {
	fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> RoundKeys<NRM1> {
		RoundKeys {
			rounds: RIArr::from_fn(|_| rng.random()),
			last: (rng.random(), rng.random()),
		}
	}
}

use crate::dual::{Dual, Q, compute_rcon};

#[instrument(level = Level::DEBUG, skip(key, sbox, dual))]
pub(crate) fn expand_nonstandard_keys<S: SBox + Copy, const NK: usize, const NRM1: usize>(
	sbox: S,
	mut key: Key<NK>,
	dual: Dual,
) -> RoundKeys<NRM1> {
	let q = Q::for_dual(dual);
	key.apply_q(q);

	const NB: usize = 4;

	let mut out = RoundKeys::new();

	debug!("initial keys");
	for i in 0..NK {
		out.set_schedule_word(i, key.word(i));
	}

	debug!("remaining keys");
	let mut rcon_idx = 1usize;
	for i in NK..(NB * (NRM1 + 2)) {
		let w_im1 = out.get_schedule_word(i - 1);
		let w_imn = out.get_schedule_word(i - NK);

		let w_i = w_imn
			^ if i.is_multiple_of(NK) {
				let rc = compute_rcon(rcon_idx, dual);
				let r = sub_word(sbox, rot_word(w_im1)) ^ Word::from_bytes([rc, 0, 0, 0]);
				rcon_idx += 1;
				r
			} else if NK > 6 && (i % NK) == 4 {
				sub_word(sbox, w_im1)
			} else {
				w_im1
			};

		out.set_schedule_word(i, w_i);
	}

	out
}
