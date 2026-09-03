/*
 * test_correctness.c -- End-to-end validation of the incremental S_q update.
 *
 * The engine maintains S_q incrementally (Troster 2008, CPC 179:30, Eq. 14).
 * After many accepted moves, any bookkeeping error would make the maintained
 * S_q drift away from the true convolution.  Here we run a batch of sweeps,
 * snapshot the maintained S, then recompute S from scratch with calcS() and
 * compare.  We also verify the reality condition h_{-q} = conj(h_q).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "membrane.h"

int main(void) {
    Config cfg = {0};
    cfg.N = 6; cfg.n = 6; cfg.p8 = 0.5; cfg.nthreads = 1;
    cfg.therm = 0; cfg.sweeps = 0; cfg.meas_every = 1; cfg.d0 = 2.6;
    cfg.seed = 999u; cfg.verbose = 0;

    Geometry geo = geometry_make(&cfg);
    int L = geo.L, LL = L * L;

    Replica rep;
    replica_alloc(&rep, &geo);
    replica_init(&rep, &geo, cfg.seed, 1u);

    for (int s = 0; s < 80; s++) replica_sweep(&rep, &geo, &cfg);

    /* snapshot the incrementally-maintained S */
    double complex *Smaintained = malloc((size_t)LL * sizeof(double complex));
    for (int i = 0; i < LL; i++) Smaintained[i] = rep.S[i];

    /* recompute from current h */
    calcS(&rep, &geo);

    double max_dS = 0.0;
    for (int i = 0; i < LL; i++) {
        double d = cabs(Smaintained[i] - rep.S[i]);
        if (d > max_dS) max_dS = d;
    }

    /* reality condition h_{-q} = conj(h_q) */
    double max_dreal = 0.0;
    for (int q1 = -cfg.N; q1 <= cfg.N; q1++)
        for (int q2 = -cfg.N; q2 <= cfg.N; q2++) {
            int i  = geo.wrap[q1 + L] * L + geo.wrap[q2 + L];
            int mi = geo.wrap[-q1 + L] * L + geo.wrap[-q2 + L];
            double d = cabs(rep.h[i] - conj(rep.h[mi]));
            if (d > max_dreal) max_dreal = d;
        }

    printf("incremental-vs-exact S : max |dS| = %.3e\n", max_dS);
    printf("reality h_{-q}=conj(h_q): max |d|  = %.3e\n", max_dreal);

    int ok = (max_dS < 1e-9) && (max_dreal < 1e-12);
    printf("%s\n", ok ? "PASS" : "FAIL");

    free(Smaintained);
    replica_free(&rep);
    geometry_free(&geo);
    return ok ? 0 : 1;
}
