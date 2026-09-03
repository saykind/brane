/*
 * pcg.h -- Minimal PCG32 pseudo-random number generator.
 *
 * A small, fast, statistically excellent RNG with independent per-stream
 * state, which is exactly what we need for replica-parallel Monte Carlo:
 * each Markov chain owns a pcg32 seeded with a distinct stream id, so the
 * chains are statistically independent and every run is fully reproducible
 * from (seed, stream).  This replaces the global, low-quality, non
 * thread-safe rand() used by the original code.
 *
 * Reference: M.E. O'Neill, "PCG: A Family of Simple Fast Space-Efficient
 * Statistically Good Algorithms for Random Number Generation" (2014),
 * https://www.pcg-random.org/  (minimal C variant, public domain / Apache-2.0).
 */
#ifndef BRANE_PCG_H
#define BRANE_PCG_H

#include <stdint.h>

typedef struct { uint64_t state; uint64_t inc; } pcg32;

static inline uint32_t pcg32_next(pcg32 *r) {
    uint64_t old = r->state;
    r->state = old * 6364136223846793005ULL + r->inc;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
}

/* Seed generator for stream `seq` (each stream is an independent sequence). */
static inline void pcg32_seed(pcg32 *r, uint64_t seed, uint64_t seq) {
    r->state = 0u;
    r->inc = (seq << 1u) | 1u;
    pcg32_next(r);
    r->state += seed;
    pcg32_next(r);
}

/* Uniform double in [0, 1). */
static inline double pcg32_double(pcg32 *r) {
    return pcg32_next(r) * (1.0 / 4294967296.0);
}

#endif /* BRANE_PCG_H */
