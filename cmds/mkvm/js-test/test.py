from forward import cipher as forward
from inv import cipher as inv

a = bytearray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
orig = a.copy()

forward(a)
if a == orig:
    raise Exception("Should not match original")
inv(a)
if a != orig:
    raise Exception("Failed to decipher")

