/*
 * test_correctness.c -- End-to-end validation of the incremental S_q update.
 *
 * The engine maintains S_q incrementally (Troster 2008, CPC 179:30, Eq. 14).
 * After many accepted moves, any bookkeeping error would make the maintained
 * S_q drift away from the true convolution.  Here we run a batch of sweeps,
 * snapshot the maintained S, then recompute S from scratch with calcS() and
 * compare.  We also verify the reality condition h_{-q} = conj(h_q).
 *
 * We run this for the plain Metropolis engine AND with over-relaxation enabled
 * (overrelax>0), since the over-relaxation move drives the same incremental S
 * machinery through a different proposal -- a bookkeeping bug there must also
 * be caught. The over-relaxation case additionally checks that the move truly
 * displaces the field (it is not a silent no-op).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "membrane.h"

/* Run `nsweeps` sweeps from a fresh replica; report the max incremental-S
 * error, max reality-condition violation, and (if `moved` != NULL) the RMS
 * field displacement from the initial configuration. */
static void check_engine(Config cfg, int nsweeps, double *max_dS,
                         double *max_dreal, double *moved) {
    Geometry geo = geometry_make(&cfg);
    int L = geo.L, LL = L * L;

    Replica rep;
    replica_alloc(&rep, &geo);
    replica_init(&rep, &geo, cfg.seed, 1u);

    double complex *h0 = NULL;
    if (moved) {
        h0 = malloc((size_t)LL * sizeof(double complex));
        for (int i = 0; i < LL; i++) h0[i] = rep.h[i];
    }

    for (int s = 0; s < nsweeps; s++) replica_sweep(&rep, &geo, &cfg);

    /* snapshot the incrementally-maintained S, then recompute from current h */
    double complex *Smaintained = malloc((size_t)LL * sizeof(double complex));
    for (int i = 0; i < LL; i++) Smaintained[i] = rep.S[i];
    calcS(&rep, &geo);

    *max_dS = 0.0;
    for (int i = 0; i < LL; i++) {
        double d = cabs(Smaintained[i] - rep.S[i]);
        if (d > *max_dS) *max_dS = d;
    }

    *max_dreal = 0.0;
    for (int q1 = -cfg.N; q1 <= cfg.N; q1++)
        for (int q2 = -cfg.N; q2 <= cfg.N; q2++) {
            int i  = geo.wrap[q1 + L] * L + geo.wrap[q2 + L];
            int mi = geo.wrap[-q1 + L] * L + geo.wrap[-q2 + L];
            double d = cabs(rep.h[i] - conj(rep.h[mi]));
            if (d > *max_dreal) *max_dreal = d;
        }

    if (moved) {
        double ss = 0.0;
        for (int i = 0; i < LL; i++) {
            double d = cabs(rep.h[i] - h0[i]);
            ss += d * d;
        }
        *moved = sqrt(ss / LL);
        free(h0);
    }

    free(Smaintained);
    replica_free(&rep);
    geometry_free(&geo);
}

int main(void) {
    Config cfg = {0};
    cfg.N = 6; cfg.n = 6; cfg.p8 = 0.5; cfg.nthreads = 1;
    cfg.therm = 0; cfg.sweeps = 0; cfg.meas_every = 1; cfg.d0 = 2.6;
    cfg.seed = 999u; cfg.verbose = 0;

    /* --- plain Metropolis --- */
    double dS_mc, dr_mc;
    check_engine(cfg, 80, &dS_mc, &dr_mc, NULL);
    printf("[metropolis]      max |dS| = %.3e   reality max |d| = %.3e\n", dS_mc, dr_mc);

    /* --- with over-relaxation --- */
    Config cor = cfg;
    cor.overrelax = 2;
    double dS_or, dr_or, moved_or;
    check_engine(cor, 80, &dS_or, &dr_or, &moved_or);
    printf("[over-relaxation] max |dS| = %.3e   reality max |d| = %.3e   rms move = %.3e\n",
           dS_or, dr_or, moved_or);

    int ok = (dS_mc < 1e-9) && (dr_mc < 1e-12)
          && (dS_or < 1e-9) && (dr_or < 1e-12)
          && (moved_or > 1e-6);   /* OR must actually displace the field */
    printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
