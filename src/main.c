/*
 * main.c -- Replica-parallel driver for the Fourier Monte Carlo membrane
 * simulation.
 *
 * Parallelization strategy: by default we run `nthreads` *independent* Markov
 * chains ("replicas"), one per core, each with its own PCG stream -- this
 * scales almost linearly and needs no per-step synchronization, and is the
 * right choice for statistics (more replicas = smaller error bars).
 *
 * For the large-N *reach* goal a single chain is the bottleneck (one sweep is
 * O(L^4) and sequential), so `inner` (it=) additionally parallelizes the
 * O(L^2) step loop across cores. Use nt=1 it=<cores> for one fast chain, or a
 * hybrid nt*it = cores. This only pays off on macOS (LLVM libomp) at large N;
 * on Linux (GCC libgomp) per-step fork/join is too costly and it regresses --
 * there, stick with replicas (it=1). Gated by lattice size (BRANE_PAR_MIN_LL).
 * See cloud/SIMCLOUD.md for measurements.
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

#ifndef GIT_SHA
#define GIT_SHA "unknown"   /* overridden by the Makefile: -DGIT_SHA=... */
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
    c.inner = 1;         /* threads per replica (intra-chain); 1 = replica-only */
    c.therm = 300;       /* thermalization sweeps (matches legacy MTH=300)   */
    c.sweeps = 2000;     /* max measurement sweeps (cap for adaptive run)   */
    c.meas_every = 1;
    c.d0 = 2.6;
    c.seed = 12345u;
    c.verbose = 0;
    c.eps = 0.005;       /* target rel. stat. error on Delta2 (0 disables)  */
    c.min_sweeps = 200;  /* don't stop before this many measurement sweeps   */
    c.block = 20;        /* sweeps between convergence checks                */
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
    printf("  nt=<int>          replicas / threads               (%d)\n", d->nthreads);
    printf("  it=<int>          threads per replica (large-N)     (%d)\n", d->inner);
    printf("  therm=<int>       thermalization sweeps per replica(%ld)\n", d->therm);
    printf("  sweeps=<int>      MAX measurement sweeps (cap)     (%ld)\n", d->sweeps);
    printf("  eps=<float>       target rel. error on Delta2; 0=off (%.3f)\n", d->eps);
    printf("  minsweeps=<int>   min sweeps before stopping       (%ld)\n", d->min_sweeps);
    printf("  block=<int>       sweeps between convergence checks(%d)\n", d->block);
    printf("  meas=<int>        measure every M sweeps           (%d)\n", d->meas_every);
    printf("  d0=<float>        base step size                   (%.2f)\n", d->d0);
    printf("  seed=<int>        base RNG seed                    (%llu)\n",
           (unsigned long long)d->seed);
    printf("  out=<path>        explicit output .dat path (overrides outdir layout)\n");
    printf("  outdir=<dir>      base dir; writes <dir>/N<N>/p<p8>/<stop>/therm..._nt..._it..._seed....dat\n");
    printf("                    <stop> = eps<eps> (adaptive) or fixed<sweeps> (eps=0)\n");
}

/* Write G(q)/G^{-1}(q) to outpath atomically (temp file + rename), so a run
 * that is killed mid-write -- e.g. by a Simcloud job timeout -- never leaves a
 * truncated file, and the previous checkpoint survives. Called both per-block
 * (checkpointing) and at the end. The header records the full run config as
 * key=value tokens (tools/analyze.py reads any k=v it finds). */
static void write_result_file(const char *outpath, const Config *cfg,
                              const Geometry *geo, const Result *res,
                              int L, long done, int converged, double wall_s) {
    char tmp[600];
    snprintf(tmp, sizeof(tmp), "%s.tmp", outpath);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", tmp); return; }
    long steps_per_sweep = (long)(2 * cfg->n + 1) * (2 * cfg->n + 1);
    fprintf(f, "# Fourier MC membrane\n");
    fprintf(f, "# N=%d L=%d n=%d p8=%.4f N8=%d Y=%.6f d0=%.4f seed=%llu\n",
            cfg->N, L, cfg->n, cfg->p8, geo->N8, geo->Y, cfg->d0,
            (unsigned long long)cfg->seed);
    fprintf(f, "# nt=%d it=%d cores=%d\n",
            cfg->nthreads, cfg->inner, cfg->nthreads * cfg->inner);
    fprintf(f, "# therm=%ld sweeps=%ld sweeps_cap=%ld min_sweeps=%ld block=%d meas_every=%d "
               "steps_per_sweep=%ld\n",
            cfg->therm, done, cfg->sweeps, cfg->min_sweeps, cfg->block,
            cfg->meas_every, steps_per_sweep);
    fprintf(f, "# eps=%.6f rel_err=%.6f converged=%d\n",
            cfg->eps, res->rel_err, converged);
    fprintf(f, "# samples=%ld accept_rate=%.4f wall_s=%.2f nu=%.6f nu_err=%.6f\n",
            res->total_meas, res->accept_rate, wall_s, res->poisson, res->poisson_err);
    fprintf(f, "# engine_sha=%s\n", GIT_SHA);
    fprintf(f, "# q1 q2 qx qy qmag G Gerr Ginv\n");
    for (int q1 = -cfg->N; q1 <= cfg->N; q1++)
        for (int q2 = -cfg->N; q2 <= cfg->N; q2++) {
            if (q1 == 0 && q2 == 0) continue;
            int i1 = geo->wrap[q1 + L], i2 = geo->wrap[q2 + L];
            double qx = geo->a * q1, qy = geo->a * q2;
            double qm = sqrt(qx * qx + qy * qy);
            double G = res->G[i1 * L + i2];
            double Ge = res->Gerr[i1 * L + i2];
            fprintf(f, "%d\t%d\t%.8f\t%.8f\t%.8f\t%.10e\t%.10e\t%.10e\n",
                    q1, q2, qx, qy, qm, G, Ge, G > 0 ? 1.0 / G : 0.0);
        }
    fclose(f);
    rename(tmp, outpath);
}

int main(int argc, char *argv[]) {
    Config cfg = default_config();
    char outpath[512] = {0};
    char outdir[400] = "data";      /* base dir; descriptive subpath appended */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(&cfg); return 0; }
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) { cfg.verbose = 1; continue; }
        if (argv[i][0] == '-') { printf("Unknown option '%s'\n", argv[i]); return 1; }
        char sbuf[512];
        if (sscanf(argv[i], "N=%d", &cfg.N) == 1) continue;
        if (sscanf(argv[i], "n=%d", &cfg.n) == 1) continue;
        if (sscanf(argv[i], "p8=%lf", &cfg.p8) == 1) continue;
        if (sscanf(argv[i], "nt=%d", &cfg.nthreads) == 1) continue;
        if (sscanf(argv[i], "it=%d", &cfg.inner) == 1) continue;
        if (sscanf(argv[i], "therm=%ld", &cfg.therm) == 1) continue;
        if (sscanf(argv[i], "sweeps=%ld", &cfg.sweeps) == 1) continue;
        if (sscanf(argv[i], "eps=%lf", &cfg.eps) == 1) continue;
        if (sscanf(argv[i], "minsweeps=%ld", &cfg.min_sweeps) == 1) continue;
        if (sscanf(argv[i], "block=%d", &cfg.block) == 1) continue;
        if (sscanf(argv[i], "meas=%d", &cfg.meas_every) == 1) continue;
        if (sscanf(argv[i], "d0=%lf", &cfg.d0) == 1) continue;
        if (sscanf(argv[i], "seed=%llu", (unsigned long long *)&cfg.seed) == 1) continue;
        if (sscanf(argv[i], "out=%511s", sbuf) == 1) { strncpy(outpath, sbuf, sizeof(outpath) - 1); continue; }
        if (sscanf(argv[i], "outdir=%399s", sbuf) == 1) { strncpy(outdir, sbuf, sizeof(outdir) - 1); continue; }
        printf("Unrecognized argument '%s'\n", argv[i]);
        return 1;
    }
    if (cfg.n <= 0 || cfg.n > cfg.N) cfg.n = cfg.N;
    if (cfg.nthreads < 1) cfg.nthreads = 1;
    if (cfg.inner < 1) cfg.inner = 1;
    if (cfg.meas_every < 1) cfg.meas_every = 1;
    omp_set_num_threads(cfg.nthreads);
    if (cfg.inner > 1) omp_set_max_active_levels(2);  /* allow nested inner teams */

    Geometry geo = geometry_make(&cfg);
    int L = geo.L;

    printf("brane: N=%d L=%d n=%d p8=%.3f  replicas=%d inner=%d (cores=%d)\n",
           cfg.N, L, cfg.n, cfg.p8, cfg.nthreads, cfg.inner, cfg.nthreads * cfg.inner);
    printf("       therm=%ld sweeps<=%ld eps=%.3f minsweeps=%ld block=%d N8=%d\n",
           cfg.therm, cfg.sweeps, cfg.eps, cfg.min_sweeps, cfg.block, geo.N8);

    /* Resolve the output path and create its directory up front, so per-block
     * checkpoints during the run can write to it. Descriptive layout keeps
     * different configs from overwriting each other:
     *   <outdir>/N<N>/p<p8>/<stop>/therm<T>_nt<NT>_it<IT>_seed<SEED>.dat
     * where <stop> = eps<eps> (adaptive) or fixed<sweeps> (fixed-length). */
    if (!outpath[0]) {
        char stop[64], fname[256];
        if (cfg.eps > 0) snprintf(stop, sizeof stop, "eps%g", cfg.eps);
        else             snprintf(stop, sizeof stop, "fixed%ld", cfg.sweeps);
        snprintf(fname, sizeof fname, "therm%ld_nt%d_it%d_seed%llu.dat",
                 cfg.therm, cfg.nthreads, cfg.inner, (unsigned long long)cfg.seed);
        snprintf(outpath, sizeof outpath, "%s/N%d/p%.2f/%s/%s",
                 outdir, cfg.N, cfg.p8, stop, fname);
    }
    {
        char dir[512];
        snprintf(dir, sizeof(dir), "%s", outpath);
        char *slash = strrchr(dir, '/');
        if (slash) {
            *slash = '\0';
            char cmd[600];
            snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", dir);
            if (system(cmd)) { /* ignore */ }
        }
    }

    Replica *reps = calloc((size_t)cfg.nthreads, sizeof(Replica));
    double t0 = omp_get_wtime();

    /* Thermalize all replicas (parallel), then measure in blocks. */
    #pragma omp parallel num_threads(cfg.nthreads)
    {
        int r = omp_get_thread_num();
        Replica *rep = &reps[r];
        replica_alloc(rep, &geo);
        replica_init(rep, &geo, cfg.seed, (uint64_t)(r + 1));
        for (long s = 0; s < cfg.therm; s++)
            replica_sweep(rep, &geo, &cfg);
    }

    /* Adaptive measurement: keep running blocks of sweeps until the relative
     * statistical error of Delta2 (estimated from the spread across the
     * independent replicas) drops below cfg.eps, bounded by [min_sweeps,
     * sweeps]. eps <= 0 disables the check (fixed cfg.sweeps sweeps). */
    long done = 0;
    double rel = -1.0;
    int converged = 0;
    int block = cfg.block > 0 ? cfg.block : 20;
    double last_ckpt = omp_get_wtime();
    const double CKPT_INTERVAL = 60.0;  /* seconds between checkpoints */

    /* Convergence trace: one row per block written to a sibling <out>.trace
     * file (flushed each block), so the sweeps/Delta2/rel_err trajectory is
     * captured on disk regardless of stdout buffering. */
    char tracepath[600];
    snprintf(tracepath, sizeof tracepath, "%s.trace", outpath);
    FILE *trace = fopen(tracepath, "w");
    if (trace) {
        fprintf(trace, "# convergence trace: N=%d p8=%.4f nt=%d it=%d therm=%ld\n",
                cfg.N, cfg.p8, cfg.nthreads, cfg.inner, cfg.therm);
        fprintf(trace, "# sweeps\tDelta2\trel_err\twall_s\n");
        fflush(trace);
    }

    while (done < cfg.sweeps) {
        long todo = block;
        if (done + todo > cfg.sweeps) todo = cfg.sweeps - done;
        #pragma omp parallel for num_threads(cfg.nthreads) schedule(static, 1)
        for (int r = 0; r < cfg.nthreads; r++) {
            Replica *rep = &reps[r];
            for (long s = 0; s < todo; s++) {
                replica_sweep(rep, &geo, &cfg);
                if ((s % cfg.meas_every) == 0) replica_measure(rep, &geo);
            }
        }
        done += todo;
        rel = delta2_rel_error(reps, cfg.nthreads, &geo);
        double d2m = delta2_mean(reps, cfg.nthreads, &geo);
        if (trace) {
            fprintf(trace, "%ld\t%.8e\t%.6f\t%.2f\n",
                    done, d2m, rel, omp_get_wtime() - t0);
            fflush(trace);
        }
        if (cfg.verbose)
            printf("  sweeps=%ld  Delta2=%.6f  Delta2 rel.err=%.4f  (target %.4f)\n",
                   done, d2m, rel, cfg.eps);
        if (cfg.eps > 0 && done >= cfg.min_sweeps && rel >= 0 && rel < cfg.eps) {
            converged = 1;
            break;
        }
        /* Checkpoint: persist current G so a killed run (e.g. job timeout)
         * still leaves usable data. Time-gated so fast small-N runs don't
         * thrash I/O. */
        if (omp_get_wtime() - last_ckpt > CKPT_INTERVAL) {
            Result cres = result_reduce(reps, cfg.nthreads, &geo);
            cres.rel_err = rel;
            write_result_file(outpath, &cfg, &geo, &cres, L, done, converged,
                              omp_get_wtime() - t0);
            result_free(&cres);
            last_ckpt = omp_get_wtime();
            if (cfg.verbose) printf("  [checkpoint @ %ld sweeps]\n", done);
        }
    }
    if (trace) fclose(trace);

    double elapsed = omp_get_wtime() - t0;
    Result res = result_reduce(reps, cfg.nthreads, &geo);
    res.sweeps_done = done;
    res.rel_err = rel;
    res.converged = converged;

    printf("\ntime = %.2f s   accept = %.1f%%   samples = %ld   sweeps/replica = %ld\n",
           elapsed, 100.0 * res.accept_rate, res.total_meas, done);
    printf("Delta2 rel.err = %.4f  %s(target %.4f)\n", res.rel_err,
           converged ? "converged " : (cfg.eps > 0 ? "NOT converged " : "fixed-length "),
           cfg.eps);
    printf("Poisson ratio nu = %.4f +/- %.4f\n", res.poisson, res.poisson_err);

    /* ---- write final Green function G(q) and inverse Green G^{-1}(q) --- */
    write_result_file(outpath, &cfg, &geo, &res, L, done, converged, elapsed);
    printf("wrote %s\n", outpath);

    for (int r = 0; r < cfg.nthreads; r++) replica_free(&reps[r]);
    free(reps);
    result_free(&res);
    geometry_free(&geo);
    return 0;
}
