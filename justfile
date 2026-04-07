run *args='-h':
	cargo run --release --bin riggle -- {{args}}

stress *args='-h':
	cargo run --release --bin stress -- {{args}}

gen *args='-h':
	cargo run --release --bin gen -- {{args}}

bench:
	cargo bench
