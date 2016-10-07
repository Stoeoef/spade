.PHONY: all
all:
	cargo build

.PHONY: update
update:
	cargo update
	make -C examples update
	rm -rf doc

.PHONY: test
test:
	cargo test

.PHONY: clean
clean:
	cargo clean
	make -C examples clean
	rm -rf doc

.PHONY: travis
travis: test
	make -C examples
