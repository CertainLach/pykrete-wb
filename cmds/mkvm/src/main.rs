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
use pykrete_wb::encoding::ExternalEncoding;
use pykrete_wb::key::{Aes128Key, Aes192Key, Aes256Key, RoundKeys};
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

#[derive(Parser)]
enum Opts {
	CreateExternalEncoding {
		#[clap(long)]
		out: String,
	},
	ExpandRoundKeys {
		t: AesTy,
		#[clap(long)]
		out: String,
		key: String,
	},
	CreateRoundKeys {
		t: AesTy,
		#[clap(long)]
		out: String,
	},
	CreateWhitebox {
		t: AesTy,
		#[clap(long)]
		key: String,

		// TODO: Provide better aliases (decode-input, encode-output) back
		#[clap(long)]
		input_encoding: Vec<String>,
		#[clap(long)]
		output_encoding: Vec<String>,

		#[clap(long)]
		inv: bool,

		#[clap(long)]
		language: Language,

		/// Pretty-print output
		#[clap(long)]
		debug: bool,
	},
}

macro_rules! match_aes_ty {
	($t:ident $(.key($key:ident))? $(.tables($tables:ident))? => $($tt:tt)*) => {
		match $t {
			AesTy::Aes128 => {
				$(type $key = Aes128Key;)?
				$(type $tables = Aes128Tables;)?
				$($tt)*
			},
			AesTy::Aes192 => {
				$(type $key = Aes192Key;)?
				$(type $tables = Aes192Tables;)?
				$($tt)*
			},
			AesTy::Aes256 => {
				$(type $key = Aes256Key;)?
				$(type $tables = Aes256Tables;)?
				$($tt)*
			}
		}
	}
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
			match_aes_ty!(t.key(AesKey) => {
				let key: RoundKeys<{AesKey::NRM1}> = rng.random();
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
		} => {
			let key = File::open(&key)?;
			match_aes_ty!(t.key(AesKey).tables(AesTables) => {
				let round_keys = RoundKeys::<{AesKey::NRM1}>::read_from(key)?;
				let mut tables = AesTables::from_unchecked_round_keys(
					&round_keys,
					inv
				);
				tables.apply_security(Security::full(), &mut rng);
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
				vmout(&tables, match language {
					Language::Js => vm::Language::Js,
					Language::Python => vm::Language::Python,
					Language::Rust => vm::Language::Rust,
					Language::Java => vm::Language::Java,
				}, debug);
			})
		}
		Opts::ExpandRoundKeys { t, out, key } => {
			let out = File::create(out)?;
			match_aes_ty!(t.key(AesKey) => {
				let key = AesKey::new(key.as_bytes().try_into().map_err(|e| anyhow::anyhow!("failed to parse key: {e}"))?);
				let rk = key.expand();
				rk.write_to(out)?;
			})
		}
	}

	Ok(())
}
