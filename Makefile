# Makefile for the brane Fourier Monte Carlo membrane simulation.
#
# OpenMP on macOS: Apple's clang needs the Homebrew `libomp` package
#   brew install libomp
# and the flags below.  On Linux with GCC, plain `-fopenmp` is enough;
# override on the command line, e.g.:  make CC=gcc-14 OMP=-fopenmp

UNAME_S := $(shell uname -s)

CC      ?= clang
CSTD     = -std=c11
WARN     = -Wall -Wextra
OPT      = -O3 -march=native -ffast-math -funroll-loops
CFLAGS   = $(CSTD) $(WARN) $(OPT) -Isrc
LDLIBS   = -lm

# --- OpenMP autodetection ------------------------------------------------
ifeq ($(UNAME_S),Darwin)
  LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(LIBOMP),)
    OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    OMP_LDLIBS = -L$(LIBOMP)/lib -lomp
  else
    $(warning libomp not found; run 'brew install libomp' for multithreading)
    OMP_CFLAGS =
    OMP_LDLIBS =
  endif
else
  OMP_CFLAGS = -fopenmp
  OMP_LDLIBS = -fopenmp
endif

CFLAGS  += $(OMP_CFLAGS)
LDLIBS  += $(OMP_LDLIBS)

SRC   = src/membrane.c
BIN   = brane
TEST  = test_correctness

.PHONY: all test clean run

all: $(BIN)

$(BIN): src/main.c $(SRC) src/membrane.h src/pcg.h
	$(CC) $(CFLAGS) src/main.c $(SRC) -o $(BIN) $(LDLIBS)

$(TEST): tests/test_correctness.c $(SRC) src/membrane.h src/pcg.h
	$(CC) $(CFLAGS) tests/test_correctness.c $(SRC) -o $(TEST) $(LDLIBS)

test: $(TEST)
	./$(TEST)

run: $(BIN)
	./$(BIN) -v

clean:
	rm -f $(BIN) $(TEST)
