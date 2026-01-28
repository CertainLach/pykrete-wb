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

use std::fs::File;

use clap::{Parser, ValueEnum};
use pykrete_wb::dual::{Dual, IRREDUCIBLE_POLYNOMIALS};
use pykrete_wb::encoding::ExternalEncoding;
use pykrete_wb::karroumi::{
	Dual4, KarroumiConfig, KarroumiConfig4, PrecomputedSBoxes, PrecomputedSBoxes4,
};
use pykrete_wb::key::{Aes128Key, Aes192Key, Aes256Key, RoundKeys};
use pykrete_wb::sbox::PrecomputedSBox;
use pykrete_wb::vm::{self, vmout};
use pykrete_wb::{Aes128Tables, Aes192Tables, Aes256Tables, Security};
use rand::{Rng, rng};

#[derive(ValueEnum, Clone, Copy, Debug)]
enum AesTy {
	Aes128,
	Aes192,
	Aes256,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum Language {
	Js,
	Python,
	Rust,
	Java,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum Karroumi {
	/// Each round uses a different AES dual, slow
	PerRound,
	/// Each row in each round uses a different AES dual, very slow
	PerRow,
}

#[derive(Parser)]
enum Opts {
	/// Create external encoding
	///
	/// Applying those will make the output incompatible with standard json, but it makes it much harder
	/// to extract the original encryption keys (given that no one has access to the external encodings).
	CreateExternalEncoding {
		/// Output file
		#[clap(long)]
		out: String,
	},
	/// Expand aes key into round keys
	///
	/// First round key is equal to the specified key, other round keys
	/// are derived from it in a predictable way.
	ExpandRoundKeys {
		/// AES standard size
		t: AesTy,
		/// Output file
		#[clap(long)]
		out: String,
		/// AES key: 16 bytes
		key: String,
	},
	/// Create round keys directly
	///
	/// Little bit secure for white-box context, as it makes it impossible to
	/// use knowledge of ExpandKeys algorithm to get the original key from one of the round keys
	CreateRoundKeys {
		/// AES standard size
		t: AesTy,
		/// Output file
		#[clap(long)]
		out: String,
	},

	/// Generate AES whitebox code
	CreateWhitebox {
		/// AES standard size
		t: AesTy,

		/// Expanded key file - can be created using expand-round-keys or create-round-keys
		#[clap(long)]
		key: String,

		/// Apply external encoding on the input, external encoding can be created using `create-external-encoding`
		/// subcommand
		///
		/// The result will not be compatible with standard AES, you would need server side support
		/// for encrypting/decrypting the input/output.
		///
		/// If prefixed with ! - the inverse of the encoding is applied
		#[clap(long)]
		input_encoding: Vec<String>,
		/// Apply external encoding on the output, external encoding can be created using `create-external-encoding`
		/// subcommand
		///
		/// The result will not be compatible with standard AES, you would need server side support
		/// for encrypting/decrypting the input/output.
		///
		/// If prefixed with ! - the inverse of the encoding is applied
		#[clap(long)]
		output_encoding: Vec<String>,

		/// If not set - AES is performing the encryption operation, otherwise - decryption
		#[clap(long)]
		inv: bool,

		/// For which language the whitebox code should be generated
		#[clap(long)]
		language: Language,

		/// Pretty-print output
		#[clap(long)]
		debug: bool,

		/// Do not reorder VM instructions
		#[clap(long)]
		no_reorder: bool,

		/// Do not use mixing bijections, reduces the attack complexity
		///
		/// Disabling both MB and L reduces code size in half
		#[clap(long)]
		no_mb: bool,
		/// Do not use mixing bijections, reduces the attack complexity
		///
		/// Disabling both MB and L reduces code size in half
		#[clap(long)]
		no_l: bool,
		/// Do not use mixing bijections, reduces the attack complexity
		///
		/// Greately reduces the code size
		#[clap(long)]
		no_internal_encodings: bool,

		/// Use Karroumi whitebox scheme instead of plain Chow
		///
		/// Greately increases VM generation time
		#[clap(long)]
		karroumi: Option<Karroumi>,

		/// Generate an AES using non-standard affine transform, MixColumns and SBox are affected
		///
		/// The result will not be compatible with standard AES, you would need server side support
		/// for encrypting/decrypting the input/output.
		///
		/// Accepts a single integer 0-240 for the given AES dual cipher, defaults to zero for standard AES
		#[clap(long, default_value = "0")]
		aes_dual: u8,
	},
}

macro_rules! match_aes_ty {
	($t:ident .key($key:ident) .tables($tables:ident), $b:block) => {
		match $t {
			AesTy::Aes128 => {
				type $key = Aes128Key;
				type $tables = Aes128Tables;
				$b
			}
			AesTy::Aes192 => {
				type $key = Aes192Key;
				type $tables = Aes192Tables;
				$b
			}
			AesTy::Aes256 => {
				type $key = Aes256Key;
				type $tables = Aes256Tables;
				$b
			}
		}
	};
}

fn main() -> anyhow::Result<()> {
	let opts = Opts::parse();
	let mut rng = rng();

	match opts {
		Opts::CreateExternalEncoding { out } => {
			let encoding: ExternalEncoding = rng.random();
			let f = File::create(out)?;
			encoding.write_to(f)?;
		}
		Opts::CreateRoundKeys { t, out } => {
			let f = File::create(out)?;
			match_aes_ty!(t.key(AesKey).tables(_T), {
				let key: RoundKeys<{ AesKey::NRM1 }> = rng.random();
				key.write_to(f)?;
			})
		}
		Opts::CreateWhitebox {
			t,
			key,
			input_encoding,
			output_encoding,
			inv,
			language,
			debug,
			no_reorder,
			aes_dual,
			no_mb,
			no_l,
			no_internal_encodings,
			karroumi,
		} => {
			let key = File::open(&key)?;
			match_aes_ty!(t.key(AesKey).tables(AesTables), {
				let round_keys = RoundKeys::<{ AesKey::NRM1 }>::read_from(key)?;
				let dual = Dual::new(
					IRREDUCIBLE_POLYNOMIALS[(aes_dual / 8) as usize],
					(aes_dual % 8) as usize,
				);
				let mut tables = match karroumi {
					Some(Karroumi::PerRound) => {
						let config = KarroumiConfig::random(&mut rng);
						let sboxes = PrecomputedSBoxes::precompute(&config);
						AesTables::from_karroumi_round_keys(
							&sboxes,
							&round_keys,
							inv,
							dual,
							&config,
						)
					}
					Some(Karroumi::PerRow) => {
						let config = KarroumiConfig4::random(&mut rng);
						let sboxes = PrecomputedSBoxes4::precompute(&config);
						AesTables::from_karroumi4_round_keys(
							&sboxes,
							&round_keys,
							inv,
							// TODO: Allow a per-row dual to be specified in cmd?
							Dual4::uniform(dual),
							&config,
						)
					}
					None => AesTables::from_nonstandard_round_keys(
						// TODO: Allow custom SBox?
						&PrecomputedSBox::aes_standard(),
						&round_keys,
						inv,
						dual,
					),
				};
				if !(no_mb && no_l && no_internal_encodings) {
					tables.apply_security(
						Security {
							mb: !no_mb,
							l: !no_l,
							internal_encodings: !no_internal_encodings,
							force_mbl: false,
							force_mbl_xor: false,
							force_xor: false,
						},
						&mut rng,
					);
				}
				for mut input in input_encoding {
					let inv = if input.starts_with("!") {
						input.remove(0);
						true
					} else {
						false
					};
					let input = ExternalEncoding::read_from(File::open(input)?)?;
					tables.input_encoding(&input, inv);
				}
				for mut output in output_encoding {
					let inv = if output.starts_with("!") {
						output.remove(0);
						true
					} else {
						false
					};
					let output = ExternalEncoding::read_from(File::open(output)?)?;
					tables.output_encoding(&output, inv);
				}
				vmout(
					&tables,
					match language {
						Language::Js => vm::Language::Js,
						Language::Python => vm::Language::Python,
						Language::Rust => vm::Language::Rust,
						Language::Java => vm::Language::Java,
					},
					debug,
					no_reorder,
				);
			})
		}
		Opts::ExpandRoundKeys { t, out, key } => {
			let out = File::create(out)?;
			match_aes_ty!(t.key(AesKey).tables(_T), {
				let key = AesKey::new(
					key.as_bytes()
						.try_into()
						.map_err(|e| anyhow::anyhow!("failed to parse key: {e}"))?,
				);
				let rk = key.expand();
				rk.write_to(out)?;
			})
		}
	}

	Ok(())
}
