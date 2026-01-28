import { cipher } from "./cipher.js";
import { cipher as forward } from "./forward.js";
import { cipher as inv } from "./inv.js";
import { cipher as encodedCipher } from "./encodedCipher.js";

import inputEncoding_ from "./input.enc" with { type: "bytes" };
import outputEncoding_ from "./output.enc" with { type: "bytes" };

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
