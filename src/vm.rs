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

use std::collections::{HashMap, HashSet};

use rand::seq::SliceRandom;
use rand::{Rng, rng};

use crate::consts::{INV_SHIFT_ROWS_TAB, SHIFT_ROWS_TAB};
use crate::{
	HighLow, Purpose, PurposeMap, RI, RowMap, SColumn, SPos, SRow, StateMap, Step, Tables, X,
};

#[derive(Clone, Copy)]
pub enum Language {
	Js,
	Python,
	Rust,
	Java,
}

struct Memory {
	id: MemoryId,
	data: Vec<u8>,
	cont_offset: usize,
}
#[derive(Default)]
struct VM {
	locals: usize,
	memory_id: usize,
	memory: Vec<Memory>,

	ops: Vec<Op>,
}
#[derive(Clone, Copy, PartialEq)]
struct MemoryId(usize);
#[derive(Clone, Copy, Hash, PartialEq, Eq, Debug)]
struct Local(usize);
impl VM {
	fn alloc_local_raw(&mut self) -> Local {
		let i = self.locals;
		self.locals += 1;
		Local(i)
	}
	fn load_state(&mut self, idx: usize) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::LoadState { idx, to });
		to
	}
	fn store_state(&mut self, from: Local, to_idx: usize) {
		self.ops.push(Op::StoreState { from, to_idx })
	}
	fn alloc_memory(&mut self, data: Vec<u8>) -> MemoryId {
		let i = self.memory_id;
		self.memory_id += 1;
		self.memory.push(Memory {
			id: MemoryId(i),
			data,
			cont_offset: 0,
		});
		MemoryId(i)
	}
	fn select(&mut self, mem: MemoryId, idx: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::Select { mem, idx, to });
		to
	}
	fn div2(&mut self, v: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::Div2 { from: v, to });
		to
	}
	fn mod2(&mut self, v: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::Mod2 { from: v, to });
		to
	}
	fn select_nib(&mut self, mem: MemoryId, idx: Local) -> Local {
		let byte_idx = self.div2(idx);
		let byte = self.select(mem, byte_idx);
		let nib_idx = self.mod2(idx);
		self.nibble_dyn(byte, nib_idx)
	}
	fn concat(&mut self, high: Local, low: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::ConcatNimbles { high, low, to });
		to
	}
	fn nibble(&mut self, from: Local, nib: HighLow) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::Nibble { from, nib, to });
		to
	}
	fn cascade1(&mut self, high: Local, low: Local, nib: HighLow) -> Local {
		let a = self.nibble(high, nib);
		let b = self.nibble(low, nib);
		self.concat(a, b)
	}
	fn xor4(&mut self, a: Local, b: Local, c: Local, d: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::Xor4 { a, b, c, d, to });
		to
	}
	fn nibble_dyn(&mut self, from: Local, nib: Local) -> Local {
		let to = self.alloc_local_raw();
		self.ops.push(Op::NibbleDyn { from, nib, to });
		to
	}
	fn finalize_memory(&mut self, no_reorder: bool) {
		if !no_reorder {
			self.memory.shuffle(&mut rng());
		}
		let mut off = 0;
		for ele in self.memory.iter_mut() {
			ele.cont_offset = off;
			off += ele.data.len();
		}
	}
}

#[derive(Clone, Copy)]
enum Op {
	LoadState {
		idx: usize,
		to: Local,
	},
	StoreState {
		from: Local,
		to_idx: usize,
	},
	Select {
		mem: MemoryId,
		idx: Local,
		to: Local,
	},
	ConcatNimbles {
		high: Local,
		low: Local,
		to: Local,
	},
	Nibble {
		from: Local,
		nib: HighLow,
		to: Local,
	},
	// to = 1 for high,
	// to = 0 for low
	NibbleDyn {
		from: Local,
		nib: Local,
		to: Local,
	},
	Div2 {
		from: Local,
		to: Local,
	},
	Mod2 {
		from: Local,
		to: Local,
	},
	Xor4 {
		a: Local,
		b: Local,
		c: Local,
		d: Local,
		to: Local,
	},
}
impl Op {
	fn into_data(self) -> OpData {
		match self {
			Op::LoadState { idx: _, to } => OpData {
				uses: [].into(),
				provides: [to].into(),
				op: self,
			},
			Op::StoreState { from, to_idx: _ } => OpData {
				uses: [from].into(),
				provides: [].into(),
				op: self,
			},
			Op::Select { mem: _, idx, to } => OpData {
				uses: [idx].into(),
				provides: [to].into(),
				op: self,
			},
			Op::ConcatNimbles { high, low, to } => OpData {
				uses: [high, low].into(),
				provides: [to].into(),
				op: self,
			},
			Op::Nibble { from, nib: _, to } => OpData {
				uses: [from].into(),
				provides: [to].into(),
				op: self,
			},
			Op::NibbleDyn { from, nib, to } => OpData {
				uses: [from, nib].into(),
				provides: [to].into(),
				op: self,
			},
			Op::Div2 { from, to } => OpData {
				uses: [from].into(),
				provides: [to].into(),
				op: self,
			},
			Op::Mod2 { from, to } => OpData {
				uses: [from].into(),
				provides: [to].into(),
				op: self,
			},
			Op::Xor4 { a, b, c, d, to } => OpData {
				uses: [a, b, c, d].into(),
				provides: [to].into(),
				op: self,
			},
		}
	}
}
struct OpData {
	uses: HashSet<Local>,
	provides: HashSet<Local>,
	op: Op,
}
#[derive(Default)]
struct LocalData {
	last_used_at: usize,
}

pub fn vmout<const NRM1: usize>(
	tables: &Tables<NRM1>,
	lang: Language,
	debug: bool,
	no_reorder: bool,
) {
	let mut vm = VM::default();
	let mut state = StateMap::from_fn(|i| vm.load_state(i.as_index()));

	let shift_rows = if tables.inv {
		&INV_SHIFT_ROWS_TAB
	} else {
		&SHIFT_ROWS_TAB
	};

	for r in RI::all::<NRM1>() {
		state = StateMap::from_fn(|pos| state[shift_rows.map(pos)]);
		for row in SRow::all() {
			for step in Step::ALL {
				if matches!(step, Step::Mbl) && !tables.uses_mbl {
					continue;
				}
				let work = match step {
					Step::Tybox => &tables.tyboxes.0[r],
					Step::Mbl => &tables.mbl.0[r],
				};
				let [aa, bb, cc, dd] = SColumn::ALL.map(|column| {
					let pos = SPos::row_column(row, column);
					RowMap::from_fn(|r| {
						let m =
							vm.alloc_memory(work.0[pos].0.into_iter().map(|c| c[r].0).collect());
						vm.select(m, state[pos])
					})
				});

				let xor = match step {
					Step::Mbl if tables.uses_xor_mbl => Some(&tables.xor_mbl),
					Step::Tybox if tables.uses_xor => Some(&tables.xor),
					_ => None,
				};
				if let Some(xor) = xor {
					let xor = xor.partial_map(r, row);

					let n01 = |vm: &mut VM, v: SRow, n: HighLow| {
						let purp = PurposeMap::from_fn(|p| {
							vm.alloc_memory(
								X::ALL
									.map(|x| {
										let (h, l) = x.as_nibs();
										xor.map(p, v, n, h, l)
									})
									.chunks(2)
									.map(|ch| {
										assert_eq!(ch.len(), 2);
										let l = ch[0];
										let h = ch[1];
										X::nibs(h, l).0
									})
									.collect(),
							)
						});
						let a = vm.cascade1(aa[v], bb[v], n);
						let a = vm.select_nib(purp[Purpose::High], a);

						let b = vm.cascade1(cc[v], dd[v], n);
						let b = vm.select_nib(purp[Purpose::Low], b);

						let o = vm.concat(a, b);
						vm.select_nib(purp[Purpose::Output], o)
					};
					let n0123 = |vm: &mut VM, v: SRow| {
						let a = n01(vm, v, HighLow::High);
						let b = n01(vm, v, HighLow::Low);
						vm.concat(a, b)
					};

					for column in SRow::ALL {
						state[SPos::row_column(row, SColumn(column.0))] = n0123(&mut vm, column);
					}
				} else {
					let n0123 = |vm: &mut VM, v: SRow| vm.xor4(aa[v], bb[v], cc[v], dd[v]);
					for column in SRow::ALL {
						state[SPos::row_column(row, SColumn(column.0))] = n0123(&mut vm, column);
					}
				}
			}
		}
	}
	state = StateMap::from_fn(|pos| state[shift_rows.map(pos)]);

	let state = StateMap::from_fn(|pos| {
		let data = vm.alloc_memory(X::ALL.map(|x| tables.tboxes_last.0[x][pos].0).to_vec());
		vm.select(data, state[pos])
	});

	for (i, l) in state.0.iter().enumerate() {
		vm.store_state(*l, i);
	}

	vm.finalize_memory(no_reorder);

	let mut opsout = vec![];
	let mut pending = vm
		.ops
		.into_iter()
		.map(|v| v.into_data())
		.collect::<Vec<_>>();
	let mut provided = HashSet::new();
	if !no_reorder {
		pending.shuffle(&mut rng());
	}
	while !pending.is_empty() {
		let mut new_pending = vec![];
		let mut added = 0;
		for ele in pending {
			if ele.uses.is_subset(&provided) {
				provided.extend(ele.provides);
				opsout.push(ele.op);
				added += 1;
			} else {
				new_pending.push(ele);
			}
		}
		assert_ne!(
			added,
			0,
			"no instructions added this iteration, {} instruction(s) remains",
			new_pending.len()
		);
		pending = new_pending;
	}
	vm.ops = opsout;

	let mut local_data = <HashMap<Local, LocalData>>::new();
	for (insn, op) in vm.ops.iter().enumerate() {
		let data = op.into_data();
		for ele in data.uses {
			local_data.entry(ele).or_default().last_used_at = insn;
		}
	}

	let mut cleanup_after_insn = <HashMap<usize, HashSet<Local>>>::new();
	for (local, data) in local_data.iter() {
		cleanup_after_insn
			.entry(data.last_used_at)
			.or_default()
			.insert(*local);
	}

	let data = vm
		.memory
		.iter()
		.flat_map(|v| &v.data)
		.copied()
		.collect::<Vec<u8>>();
	let mut reg = RegAlloc::new(lang, no_reorder);

	match lang {
		Language::Js => {
			print!("const M='");
			for ele in data.chunks(2) {
				assert!(ele.len() == 2, "expected odd memory size");
				let a = ele[0];
				let b = ele[1];
				print!("\\u{a:0>2x}{b:0>2x}");
			}
			print!("'.split('').flatMap(v=>{{let x=v.charCodeAt(0);return[x>>8,x&255]}});");
			print!("export function cipher(d){{");
		}
		Language::Python => {
			print!("M=b'");
			for b in data.iter() {
				print!("\\x{b:0>2x}");
			}
			println!("'");
			println!("def cipher(d):")
		}
		Language::Rust => {
			print!("const M:[u8;{}]=*b\"", data.len());
			for b in data.iter() {
				print!("\\x{b:0>2x}");
			}
			print!("\";");
			print!("pub fn cipher(d:&mut[u8]){{")
		}
		Language::Java => {
			print!(
				"public class A{{private static byte[] M=java.util.HexFormat.of().parseHex(new StringBuilder(\""
			);
			for (i, chunk) in data.chunks(64000 / 2).enumerate() {
				if i != 0 {
					print!("\")");
					if debug {
						println!();
					}
					print!(".append(\"")
				}
				for b in chunk.iter() {
					print!("{b:0>2x}");
				}
			}
			print!("\").toString());");
			if debug {
				println!();
			}
			print!("public static void _c(byte[]d){{");
		}
	}

	let mut splits = 0;
	for (insn, op) in vm.ops.iter().enumerate() {
		if matches!(lang, Language::Java) && insn % 500 == 0 && insn != 0 {
			let (pass, params) = reg.split_function();
			print!("_c{splits}(d");
			for ele in pass {
				print!(",{ele}")
			}
			print!(");}}");
			if debug {
				println!();
			}
			print!("private static void _c{splits}(byte[]d");
			for ele in params {
				print!(",{ele}")
			}
			print!("){{");
			splits += 1
		}

		if matches!(lang, Language::Python) {
			print!("\t");
		}
		match op {
			Op::Select { mem, idx, to } => {
				let v = format!(
					"M[{}+{}]",
					vm.memory
						.iter()
						.find(|m| m.id == *mem)
						.expect("memdata")
						.cont_offset,
					reg.used_index(*idx),
				);
				print!("{}", reg.alloc(*to, v))
			}
			Op::ConcatNimbles { high, low, to } => {
				let v = format!("({}<<4)|{}", reg.used(*high), reg.used(*low));
				print!("{}", reg.alloc(*to, v));
			}
			Op::Nibble { from, nib, to } => {
				let v = if matches!(nib, HighLow::High) {
					format!("{}>>4", reg.used(*from))
				} else {
					format!("{}&15", reg.used(*from))
				};
				print!("{}", reg.alloc(*to, v))
			}
			Op::NibbleDyn { from, nib, to } => {
				let v = format!("({}>>(4*{}))&15", reg.used(*from), reg.used(*nib));
				print!("{}", reg.alloc(*to, v))
			}
			Op::Div2 { from, to } => {
				let v = format!("{}>>1", reg.used(*from));
				print!("{}", reg.alloc(*to, v))
			}
			Op::Mod2 { from, to } => {
				let v = format!("{}&1", reg.used(*from));
				print!("{}", reg.alloc(*to, v))
			}
			Op::LoadState { idx, to } => {
				print!("{}", reg.alloc(*to, format!("d[{idx}]")));
			}
			Op::StoreState { from, to_idx } => {
				if matches!(lang, Language::Java) {
					print!("d[{to_idx}]=(byte) {}", reg.used(*from))
				} else {
					print!("d[{to_idx}]={}", reg.used(*from))
				}
			}
			Op::Xor4 { a, b, c, d, to } => {
				let v = format!(
					"{}^{}^{}^{}",
					reg.used(*a),
					reg.used(*b),
					reg.used(*c),
					reg.used(*d)
				);
				print!("{}", reg.alloc(*to, v))
			}
		}
		for ele in cleanup_after_insn
			.get(&insn)
			.unwrap_or(&HashSet::new())
			.iter()
		{
			reg.free(*ele);
		}
		if matches!(lang, Language::Python) {
			println!()
		} else {
			print!(";");
			if debug {
				println!();
			}
		}
	}
	if matches!(lang, Language::Js | Language::Rust) {
		print!("}}");
	} else if matches!(lang, Language::Java) {
		print!("}}}}");
	}
}

struct RegAlloc {
	allocated: HashMap<Local, String>,
	free_list: Vec<String>,
	last_allocated: usize,
	lang: Language,
	no_reorder: bool,
}
impl RegAlloc {
	fn new(lang: Language, no_reorder: bool) -> Self {
		Self {
			allocated: HashMap::new(),
			free_list: vec![],
			last_allocated: 0,
			lang,
			no_reorder,
		}
	}

	fn split_function(&mut self) -> (Vec<String>, Vec<String>) {
		assert!(matches!(self.lang, Language::Java));
		let mut pass = vec![];
		let mut params = vec![];

		self.free_list = vec![];

		for (_, v) in &self.allocated {
			pass.push(v.clone());
			params.push(format!("int {v}"));
		}

		(pass, params)
	}

	fn alloc_name(&mut self, l: Local) -> (bool, String) {
		if self.allocated.contains_key(&l) {
			panic!("local is already allocated");
		}
		if self.free_list.is_empty() {
			let i = self.last_allocated;
			self.last_allocated += 1;
			let name = format!("r{i}");
			self.allocated.insert(l, name.clone());
			(true, name)
		} else {
			let i = if self.no_reorder {
				0
			} else {
				rng().random_range(0..self.free_list.len())
			};
			let i = self.free_list.remove(i);
			self.allocated.insert(l, i.clone());
			(false, i)
		}
	}
	fn alloc(&mut self, l: Local, v: String) -> String {
		let (new, name) = self.alloc_name(l);
		match self.lang {
			Language::Js => {
				if new {
					format!("let {name}={v}")
				} else {
					format!("{name}={v}")
				}
			}
			Language::Rust => {
				if new {
					format!("let mut {name}:u8={v}")
				} else {
					format!("{name}={v}")
				}
			}
			Language::Java => {
				if new {
					format!("int {name}={v}")
				} else {
					format!("{name}={v}")
				}
			}
			Language::Python => {
				format!("{name}={v}")
			}
		}
	}
	fn free(&mut self, l: Local) {
		let n = self.allocated.remove(&l).expect("is allocated before free");
		self.free_list.push(n);
	}

	fn used(&mut self, l: Local) -> String {
		self.allocated
			.get(&l)
			.expect("is allocated before used")
			.clone()
	}
	fn used_index(&mut self, l: Local) -> String {
		if matches!(self.lang, Language::Rust) {
			format!("{} as usize", self.used(l))
		} else {
			self.used(l)
		}
	}
}
