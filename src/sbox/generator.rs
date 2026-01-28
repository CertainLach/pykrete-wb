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

use core::array::from_fn;

use rand::prelude::*;

use crate::{X, XArr};

pub fn random_sbox<R: Rng>(rng: &mut R) -> XArr<X> {
	let mut sbox: [X; 256] = from_fn(|i| X(i as u8));
	sbox.shuffle(rng);
	let mut score = evaluate_sbox(&sbox);

	for _ in 0..1000 {
		let mut candidate = sbox;

		let i = rng.random_range(0..256);
		let j = rng.random_range(0..256);
		candidate.swap(i, j);
		let candidate_score = evaluate_sbox(&candidate);

		if candidate_score > score {
			sbox = candidate;
			score = candidate_score;
		}
	}

	XArr(sbox)
}

fn evaluate_sbox(sbox: &[X; 256]) -> i32 {
	let mut score = 0;

	for (i, v) in sbox.iter().enumerate() {
		if i == v.0 as usize {
			score -= 100;
		}
	}
	for (i, v) in sbox.iter().enumerate() {
		if i == (!v.0) as usize {
			score -= 50;
		}
	}
	score += nonlinearity_score(sbox);
	score += balance_score(sbox);
	score
}

#[allow(clippy::needless_range_loop)]
fn nonlinearity_score(sbox: &[X; 256]) -> i32 {
	let mut score = 0;
	for bit in 0..8 {
		let mask = 1u8 << bit;
		let mut linear_combinations = [0; 256];

		// Walsh-hadamard
		for a in 1..256 {
			let mut correlation = 0;
			for x in 0..256 {
				let fx = ((sbox[x].0 & mask) >> bit) as i32;
				let ax = dot_product(a as u8, x as u8) as i32;
				correlation += (-1_i32).pow((fx ^ ax) as u32);
			}
			linear_combinations[a] = correlation.abs();
		}

		// Higher nonlinearity => lower max correlation
		let max_correlation = linear_combinations.iter().max().unwrap_or(&256);
		score += 256 - max_correlation;
	}

	score / 8
}

fn balance_score(sbox: &[X; 256]) -> i32 {
	let mut score = 0;

	for bit in 0..8 {
		let mask = 1u8 << bit;
		let ones = sbox
			.iter()
			.map(|&x| ((x.0 & mask) >> bit) as i32)
			.sum::<i32>();
		let balance = (128 - (ones - 128).abs()) * 2;
		score += balance;
	}

	score
}

fn dot_product(a: u8, x: u8) -> u8 {
	(a & x).count_ones() as u8 & 1
}
