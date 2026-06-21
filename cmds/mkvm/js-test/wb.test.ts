import { cipher } from "./cipher.js";
import { cipher as forward } from "./forward.js";
import { cipher as inv } from "./inv.js";
import { cipher as encodedCipher } from "./encodedCipher.js";
import { cipher as linearCipher } from "./linearCipher.js";
import { cipher as linearCipherNoinv } from "./linearCipherNoinv.js";
import { cipher as fullCipher } from "./fullCipher.js";

import inputEncoding_ from "./input.enc" with { type: "bytes" };
import outputEncoding_ from "./output.enc" with { type: "bytes" };
import linInMatrix_ from "./in.lin" with { type: "bytes" };
import linOutMatrix_ from "./out.lin" with { type: "bytes" };

import { assertEquals, assertNotEquals } from "jsr:@std/assert";

const TWO_ONE_NINE_TWO_TEST_VECTOR = "Two One Nine Two";
const TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR = [
  41,
  195,
  80,
  95,
  87,
  20,
  32,
  246,
  64,
  34,
  153,
  179,
  26,
  2,
  215,
  58,
];

function stringToBytes(v: string): number[] {
  return v.split("").map((c) => c.charCodeAt(0));
}

Deno.test(function kungFuTestVector() {
  const o = stringToBytes(TWO_ONE_NINE_TWO_TEST_VECTOR);
  cipher(o);
  assertEquals(o, TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
});

Deno.test(function encodedBidir() {
  const a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
  const orig = [...a];

  forward(a);
  assertNotEquals(a, orig);
  inv(a);
  assertEquals(a, orig);
});

function parseEncoding(data: Uint8Array): Uint8Array[] {
  const out = [];
  assertEquals(data.length, 16 * 256);
  for (let i = 0; i < 16; i++) {
    out.push(data.slice(i * 256, (i + 1) * 256));
  }
  return out;
}
function applyEncoding(state: number[], encoding: Uint8Array[], inv: boolean) {
  if (inv) {
    for (let i = 0; i < 16; i++) {
      state[i] = encoding[i].indexOf(state[i]);
    }
  } else {
    for (let i = 0; i < 16; i++) {
      state[i] = encoding[i][state[i]];
    }
  }
}

Deno.test(function manualEncoding() {
  const inputEncoding = parseEncoding(inputEncoding_);
  const outputEncoding = parseEncoding(outputEncoding_);

  const o = stringToBytes(TWO_ONE_NINE_TWO_TEST_VECTOR);
  applyEncoding(o, inputEncoding, false);
  encodedCipher(o);
  applyEncoding(o, outputEncoding, true);
  assertEquals(o, TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
});

function parseMatrix(data: Uint8Array): number[][] {
  assertEquals(data.length, 128 * 16);
  const M: number[][] = [];
  for (let i = 0; i < 128; i++) {
    const row: number[] = [];
    for (let j = 0; j < 128; j++) {
      row.push((data[i * 16 + (j >> 3)] >> (7 - (j & 7))) & 1);
    }
    M.push(row);
  }
  return M;
}
function applyLinear(state: number[], M: number[][]) {
  const v = new Array(128).fill(0);
  for (let idx = 0; idx < 16; idx++) {
    for (let bit = 0; bit < 8; bit++) {
      v[idx * 8 + bit] = (state[idx] >> (7 - bit)) & 1;
    }
  }
  for (let i = 0; i < 16; i++) state[i] = 0;
  for (let i = 0; i < 128; i++) {
    let s = 0;
    for (let j = 0; j < 128; j++) s ^= M[i][j] & v[j];
    if (s) state[i >> 3] |= 1 << (7 - (i & 7));
  }
}

function invertMatrix(M: number[][]): number[][] {
  const n = 128;
  const a = M.map((row, i) => {
    const id = new Array(n).fill(0);
    id[i] = 1;
    return row.concat(id);
  });
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let r = col; r < n; r++) {
      if (a[r][col] === 1) {
        pivot = r;
        break;
      }
    }
    if (pivot === -1) throw new Error("singular matrix");
    if (pivot !== col) [a[pivot], a[col]] = [a[col], a[pivot]];
    for (let r = 0; r < n; r++) {
      if (r !== col && a[r][col] === 1) {
        for (let k = col; k < 2 * n; k++) a[r][k] ^= a[col][k];
      }
    }
  }
  return a.map((row) => row.slice(n));
}

Deno.test(function linearEncoding() {
  const linIn = parseMatrix(linInMatrix_);
  const linOut = parseMatrix(linOutMatrix_);

  const o = stringToBytes(TWO_ONE_NINE_TWO_TEST_VECTOR);
  applyLinear(o, linIn);
  linearCipher(o);
  applyLinear(o, linOut);
  assertEquals(o, TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
});

Deno.test(function linearEncodingNoInv() {
  const linInInv = invertMatrix(parseMatrix(linInMatrix_));
  const linOutInv = invertMatrix(parseMatrix(linOutMatrix_));

  const o = stringToBytes(TWO_ONE_NINE_TWO_TEST_VECTOR);
  applyLinear(o, linInInv);
  linearCipherNoinv(o);
  applyLinear(o, linOutInv);
  assertEquals(o, TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
});

Deno.test(function fullEncoding() {
  const inputEncoding = parseEncoding(inputEncoding_);
  const outputEncoding = parseEncoding(outputEncoding_);
  const linIn = parseMatrix(linInMatrix_);
  const linOut = parseMatrix(linOutMatrix_);

  const o = stringToBytes(TWO_ONE_NINE_TWO_TEST_VECTOR);
  // input: nonlinear, then linear (outermost)
  applyEncoding(o, inputEncoding, false);
  applyLinear(o, linIn);
  fullCipher(o);
  // output: undo linear, then nonlinear
  applyLinear(o, linOut);
  applyEncoding(o, outputEncoding, true);
  assertEquals(o, TWO_ONE_NINE_TWO_AES128_KUNG_FU_TEST_VECTOR);
});
