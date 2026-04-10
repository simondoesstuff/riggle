run *args='-h':
	cargo run --release --bin riggle -- {{args}}

stress *args='-h':
	cargo run --release --bin stress -- {{args}}

gen *args='-h':
	cargo run --release --bin gen -- {{args}}

@sort *args='-h':
	if [ ! -f ./target/release/bedsort ]; then cargo build --release --bin bedsort; fi
	./target/release/bedsort {{args}}

bench:
	cargo bench

test:
	cargo test
