# T-EVOLVE reference implementation.
# Linux/macOS: make        Windows (MSYS2 / MinGW-w64): mingw32-make, or scripts\build_windows.bat
CC        ?= gcc
OPT       ?= -O3 -march=native
GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
CFLAGS    += $(OPT) -std=gnu11 -Wall -Wextra -Wno-unused-function -Iinclude \
             '-DGIT_COMMIT="$(GIT_COMMIT)"' '-DBUILD_FLAGS="$(CC) $(OPT)"'
LDLIBS    += -lm

CORE    = src/keccak.c src/ring.c src/sample.c
PROTO   = src/params.c src/bdlop.c src/shamir.c src/codec.c src/ballot.c src/agg.c src/tally.c src/evolve.c
SRC     = $(CORE) $(PROTO)

BIN     = build/test_core build/test_protocol build/exp_ballot_rej build/bench_protocol build/param_report

all: $(BIN)

build:
	mkdir -p build

build/test_core: tests/test_core.c $(CORE) | build
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%: tests/%.c $(SRC) | build
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%: experiments/%.c $(SRC) | build
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%: bench/%.c $(SRC) | build
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%: tools/%.c $(SRC) | build
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test: build/test_core build/test_protocol
	./build/test_core
	./build/test_protocol -q

clean:
	rm -rf build

.PHONY: all test clean
