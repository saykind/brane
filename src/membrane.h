/*
 * membrane.h -- Public API and data structures for the Fourier Monte Carlo
 * simulation of a 2D crystalline membrane (graphene) in the flat phase.
 *
 * Physics (dimensionless free energy, thesis eq. in chapters/mcmc.tex):
 *
 *     F = sum_{q != 0} { 1/2 (2 pi q / L)^4 |h_q|^2
 *                        + (2 pi / 3) (p8^2 / L^2) |S_q|^2 }
 *
 *     S_q = sum_{k != 0} [ (2 pi/L) k x qhat ]^2 (h_k . conj(h_{k+q}))
 *
 * with momenta represented by their lattice sine values.  The height field
 * h_q lives on an L x L grid (L = 2N+1) and obeys the reality condition
 * h_{-q} = conj(h_q).
 *
 * All arrays are stored as contiguous, flattened, row-major L*L blocks
 * indexed by IDX(q1, q2) with q1, q2 in [0, L).  Negative momenta wrap
 * around (q -> q + L), matching a discrete Brillouin zone.
 */
#ifndef BRANE_MEMBRANE_H
#define BRANE_MEMBRANE_H

#include <complex.h>
#include <stdint.h>
#include "pcg.h"

/* ---- simulation parameters ------------------------------------------- */
typedef struct {
    int      N;          /* half linear size; lattice is L = 2N+1 across   */
    int      n;          /* half size of move-selection zone (l = 2n+1)    */
    double   p8;         /* interaction strength (0 < p8 < pi)             */
    int      nthreads;   /* number of independent replicas / OpenMP threads*/
    long     therm;      /* thermalization sweeps per replica              */
    long     sweeps;     /* measurement sweeps per replica                 */
    int      meas_every; /* take a measurement every this many sweeps      */
    double   d0;         /* base Metropolis step size                      */
    uint64_t seed;       /* base RNG seed (reproducibility)                */
    int      verbose;    /* progress reporting                             */
} Config;

/* Derived, read-only geometry shared by all replicas. */
typedef struct {
    int     N, L;
    double  a;           /* lattice spacing in q-space: 2 pi / L           */
    double  Y;           /* interaction coefficient (2 pi / 3) p8^2        */
    int     N8;          /* Poisson-ratio window half-width                */
    double *SN;          /* [L]   sin(a * q)                               */
    double *Q;           /* [L*L] lattice q^2 = 4(sin^2 + sin^2)           */
    double *invQ;        /* [L*L] 1/Q (0 at the origin)                    */
    double *dstep;       /* [L*L] per-mode Metropolis step size d(k)       */
    int    *wrap;        /* [3L]  modulo table, wrap[x + L] = x mod L      */
} Geometry;

/* Per-replica state and accumulators. */
typedef struct {
    int             L;
    double complex *h;   /* [L*L] height modes                             */
    double complex *S;   /* [L*L] nonlinear array                          */
    double complex *dS;  /* [L*L] increment scratch                        */
    double         *g;   /* [L*L] accumulator: sum of |h_q|^2 over samples */
    long            nmeas;
    /* Poisson-ratio accumulators over the small-q window                  */
    double sum_Kx, sum_Ky, sum_KxKx, sum_KxKy;
    /* Metropolis bookkeeping                                              */
    long   proposed, accepted;
    pcg32  rng;
} Replica;

/* ---- API ------------------------------------------------------------- */
Geometry geometry_make(const Config *cfg);
void     geometry_free(Geometry *geo);

void     replica_alloc(Replica *rep, const Geometry *geo);
void     replica_free(Replica *rep);
void     replica_init(Replica *rep, const Geometry *geo, uint64_t seed, uint64_t stream);

/* One full sweep = l*l attempted single-mode Metropolis updates.         */
void     replica_sweep(Replica *rep, const Geometry *geo, const Config *cfg);
/* Accumulate observables (Green function, Poisson correlators).          */
void     replica_measure(Replica *rep, const Geometry *geo);

/* Recompute S_q from scratch (O(L^4)); used after init and as a self-check*/
void     calcS(Replica *rep, const Geometry *geo);

/* Combine per-replica accumulators into global averages (called serially).*/
typedef struct {
    int     N, L, N8;
    double  a, p8;
    long    total_meas;
    double *G;           /* [L*L] Green function <|h_q|^2>                 */
    double  poisson;     /* Poisson ratio nu                              */
    double  accept_rate;
} Result;

Result   result_reduce(const Replica *reps, int nreps, const Geometry *geo);
void     result_free(Result *res);

#endif /* BRANE_MEMBRANE_H */
