/*
 * main.c -- Replica-parallel driver for the Fourier Monte Carlo membrane
 * simulation.
 *
 * Parallelization strategy: instead of parallelizing the O(L^2) inner loop
 * of a single Markov chain (the original approach, which drowned in OpenMP
 * fork/join overhead at small L), we run `nthreads` *independent* Markov
 * chains ("replicas"), one per core, each with its own PCG stream.  Every
 * chain thermalizes and samples on its own; the Green function and Poisson
 * correlators are averaged across chains at the end.  This scales almost
 * linearly and needs no per-step synchronization.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "membrane.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static Config default_config(void) {
    Config c;
    c.N = 36;
    c.n = 0;              /* 0 => set to N later */
    c.p8 = 0.4;          /* ~ graphene at T=300 K (see chapters/overview.tex) */
    /* Default to 12 replicas: on the M4 Max (12 P + 4 E cores) throughput
     * peaks near the performance-core count; the E-cores gate the end-of-run
     * barrier and add nothing (see tools/scaling.py). Capped by availability. */
    int mx = omp_get_max_threads();
    c.nthreads = mx < 12 ? mx : 12;
    c.therm = 80;
    c.sweeps = 80;
    c.meas_every = 1;
    c.d0 = 2.6;
    c.seed = 12345u;
    c.verbose = 0;
    return c;
}

static void usage(const Config *d) {
    printf("brane -- Fourier Monte Carlo of a 2D crystalline membrane\n\n");
    printf("Usage: brane [options]\n\n");
    printf("  -h, --help        show this message\n");
    printf("  -v, --verbose     progress reporting\n");
    printf("  N=<int>           half lattice size, L=2N+1        (%d)\n", d->N);
    printf("  n=<int>           half move-zone size, l=2n+1      (=N)\n");
    printf("  p8=<float>        interaction strength, 0<p8<pi    (%.2f)\n", d->p8);
    printf("  nt=<int>          replicas / threads               (%d)\n", d->nthreads);    printf("  therm=<int>       thermalization sweeps per replica(%ld)\n", d->therm);
    printf("  sweeps=<int>      measurement sweeps per replica   (%ld)\n", d->sweeps);
    printf("  meas=<int>        measure every M sweeps           (%d)\n", d->meas_every);
    printf("  d0=<float>        base step size                   (%.2f)\n", d->d0);
    printf("  seed=<int>        base RNG seed                    (%llu)\n",
           (unsigned long long)d->seed);
    printf("  out=<path>        output .dat file (default data/N=<N>.dat)\n");
}

int main(int argc, char *argv[]) {
    Config cfg = default_config();
    char outpath[512] = {0};

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(&cfg); return 0; }
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) { cfg.verbose = 1; continue; }
        if (argv[i][0] == '-') { printf("Unknown option '%s'\n", argv[i]); return 1; }
        char sbuf[512];
        if (sscanf(argv[i], "N=%d", &cfg.N) == 1) continue;
        if (sscanf(argv[i], "n=%d", &cfg.n) == 1) continue;
        if (sscanf(argv[i], "p8=%lf", &cfg.p8) == 1) continue;
        if (sscanf(argv[i], "nt=%d", &cfg.nthreads) == 1) continue;
        if (sscanf(argv[i], "therm=%ld", &cfg.therm) == 1) continue;
        if (sscanf(argv[i], "sweeps=%ld", &cfg.sweeps) == 1) continue;
        if (sscanf(argv[i], "meas=%d", &cfg.meas_every) == 1) continue;
        if (sscanf(argv[i], "d0=%lf", &cfg.d0) == 1) continue;
        if (sscanf(argv[i], "seed=%llu", (unsigned long long *)&cfg.seed) == 1) continue;
        if (sscanf(argv[i], "out=%511s", sbuf) == 1) { strncpy(outpath, sbuf, sizeof(outpath) - 1); continue; }
        printf("Unrecognized argument '%s'\n", argv[i]);
        return 1;
    }
    if (cfg.n <= 0 || cfg.n > cfg.N) cfg.n = cfg.N;
    if (cfg.nthreads < 1) cfg.nthreads = 1;
    if (cfg.meas_every < 1) cfg.meas_every = 1;
    omp_set_num_threads(cfg.nthreads);

    Geometry geo = geometry_make(&cfg);
    int L = geo.L;

    printf("brane: N=%d L=%d n=%d p8=%.3f  replicas=%d\n",
           cfg.N, L, cfg.n, cfg.p8, cfg.nthreads);
    printf("       therm=%ld sweeps=%ld meas_every=%d N8=%d\n",
           cfg.therm, cfg.sweeps, cfg.meas_every, geo.N8);

    Replica *reps = calloc((size_t)cfg.nthreads, sizeof(Replica));
    double t0 = omp_get_wtime();

    #pragma omp parallel num_threads(cfg.nthreads)
    {
        int r = omp_get_thread_num();
        Replica *rep = &reps[r];
        replica_alloc(rep, &geo);
        replica_init(rep, &geo, cfg.seed, (uint64_t)(r + 1));

        for (long s = 0; s < cfg.therm; s++)
            replica_sweep(rep, &geo, &cfg);
        for (long s = 0; s < cfg.sweeps; s++) {
            replica_sweep(rep, &geo, &cfg);
            if ((s % cfg.meas_every) == 0) replica_measure(rep, &geo);
        }
        if (cfg.verbose)
            #pragma omp critical
            printf("  replica %2d: acc=%.1f%% meas=%ld\n",
                   r, 100.0 * rep->accepted / (rep->proposed ? rep->proposed : 1),
                   rep->nmeas);
    }

    double elapsed = omp_get_wtime() - t0;
    Result res = result_reduce(reps, cfg.nthreads, &geo);

    printf("\ntime = %.2f s   accept = %.1f%%   samples = %ld\n",
           elapsed, 100.0 * res.accept_rate, res.total_meas);
    printf("Poisson ratio nu = %.4f\n", res.poisson);

    /* ---- write Green function G(q) and inverse Green G^{-1}(q) --------- */
    if (!outpath[0]) {
        system("mkdir -p data");
        snprintf(outpath, sizeof(outpath), "data/N=%d.dat", cfg.N);
    }
    FILE *f = fopen(outpath, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", outpath); return 1; }
    fprintf(f, "# Fourier MC membrane   N=%d L=%d p8=%.4f samples=%ld nu=%.6f\n",
            cfg.N, L, cfg.p8, res.total_meas, res.poisson);
    fprintf(f, "# q1 q2 qx qy qmag G Ginv\n");
    for (int q1 = -cfg.N; q1 <= cfg.N; q1++)
        for (int q2 = -cfg.N; q2 <= cfg.N; q2++) {
            if (q1 == 0 && q2 == 0) continue;
            int i1 = geo.wrap[q1 + L], i2 = geo.wrap[q2 + L];
            double qx = geo.a * q1, qy = geo.a * q2;
            double qm = sqrt(qx * qx + qy * qy);
            double G = res.G[i1 * L + i2];
            fprintf(f, "%d\t%d\t%.8f\t%.8f\t%.8f\t%.10e\t%.10e\n",
                    q1, q2, qx, qy, qm, G, G > 0 ? 1.0 / G : 0.0);
        }
    fclose(f);
    printf("wrote %s\n", outpath);

    for (int r = 0; r < cfg.nthreads; r++) replica_free(&reps[r]);
    free(reps);
    result_free(&res);
    geometry_free(&geo);
    return 0;
}
