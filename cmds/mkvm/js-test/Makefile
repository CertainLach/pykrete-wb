test: test-deps
	deno test --unstable-raw-imports

clean:
	rm -f *.enc *.lin *.key *.js inv.py forward.py

%.enc:
	cargo run create-external-encoding --out $@

%.lin:
	cargo run create-linear-encoding --out $@

data.key:
	cargo run create-round-keys aes256 --out $@
forward.js: data.key input.enc output.enc
	cargo run --release create-whitebox aes256 --key data.key --input-encoding !input.enc --output-encoding output.enc --language js > $@
inv.js: data.key input.enc output.enc
	cargo run --release create-whitebox aes256 --inv --key data.key --input-encoding !output.enc --output-encoding input.enc --language js > $@
forward.py: data.key input.enc output.enc
	cargo run --release create-whitebox aes256 --key data.key --input-encoding !input.enc --output-encoding output.enc --language python > $@
inv.py: data.key input.enc output.enc
	cargo run --release create-whitebox aes256 --inv --key data.key --input-encoding !output.enc --output-encoding input.enc --language python > $@

kung_fu.key:
	cargo run expand-round-keys aes128 "Thats my Kung Fu" --out $@
cipher.js: kung_fu.key
	cargo run --release create-whitebox aes128 --key kung_fu.key --language js > $@
encodedCipher.js: kung_fu.key input.enc output.enc
	cargo run --release create-whitebox aes128 --key kung_fu.key --input-encoding !input.enc --output-encoding output.enc --language js > $@
linearCipher.js: kung_fu.key in.lin out.lin
	cargo run --release create-whitebox aes128 --key kung_fu.key --input-linear-encoding !in.lin --output-linear-encoding !out.lin --language js > $@
linearCipherNoinv.js: kung_fu.key in.lin out.lin
	cargo run --release create-whitebox aes128 --key kung_fu.key --input-linear-encoding in.lin --output-linear-encoding out.lin --language js > $@
fullCipher.js: kung_fu.key input.enc output.enc in.lin out.lin
	cargo run --release create-whitebox aes128 --key kung_fu.key --input-encoding !input.enc --input-linear-encoding !in.lin --output-encoding output.enc --output-linear-encoding !out.lin --language js > $@

test-deps: forward.js inv.js cipher.js encodedCipher.js linearCipher.js linearCipherNoinv.js fullCipher.js inv.py forward.py
