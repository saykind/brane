/*
 * membrane.c -- Core Fourier Monte Carlo engine.
 *
 * The expensive part of the model is the anharmonic term sum_q |S_q|^2 with
 * S_q a convolution over the whole lattice.  Recomputing it from scratch
 * every move costs O(L^4).  Following Troster (Monte Carlo simulation in
 * Fourier space, 2008) we instead store S_q and, for a single-mode move
 * h_{k} -> h_{k} + d*z (and its conjugate partner h_{-k} += d*conj(z)),
 * compute only the increment dS_q analytically in O(L^2).  The energy change
 * then follows from
 *
 *     dE_anharm = -(Y/L^2) sum_q Re[ (2 S_q + dS_q) conj(dS_q) ]
 *     dE_harm   = Q_k^2 * Re[ (2 h_k + d*z) conj(d*z) ]
 *
 * (signs chosen so the Metropolis weight is exp(w), w = -dE).  The dS_q
 * formula below is a faithful port of the reference implementation, which
 * was validated term-by-term against a brute-force convolution in
 * legacy/simexact.c (agreement to 5e-11).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "membrane.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IDX(q1, q2, L) ((q1) * (L) + (q2))

/* ---- geometry -------------------------------------------------------- */
Geometry geometry_make(const Config *cfg) {
    Geometry geo;
    int N = cfg->N, L = 2 * N + 1;
    double a = 2.0 * M_PI / L;
    geo.N = N; geo.L = L; geo.a = a;
    geo.Y = (2.0 * M_PI / 3.0) * cfg->p8 * cfg->p8;
    geo.N8 = (int)(cfg->p8 / M_PI * N);

    geo.SN = malloc((size_t)L * sizeof(double));
    geo.Q = malloc((size_t)L * L * sizeof(double));
    geo.invQ = malloc((size_t)L * L * sizeof(double));
    geo.dstep = malloc((size_t)L * L * sizeof(double));
    geo.wrap = malloc((size_t)(3 * L) * sizeof(int));

    for (int i = 0; i < 3 * L; i++) geo.wrap[i] = ((i - L) % L + L) % L;

    for (int q1 = -N; q1 <= N; q1++) {
        int i1 = (q1 + L) % L;
        geo.SN[i1] = sin(a * q1);
        for (int q2 = -N; q2 <= N; q2++) {
            int i2 = (q2 + L) % L;
            double q2sq = 4.0 * (sin(a * q1 / 2) * sin(a * q1 / 2) +
                                 sin(a * q2 / 2) * sin(a * q2 / 2));
            int idx = IDX(i1, i2, L);
            geo.Q[idx] = q2sq;
            geo.invQ[idx] = (q2sq != 0.0) ? 1.0 / q2sq : 0.0;
            /* momentum-dependent step d(k) ~ 1/q^2, softened by interaction */
            geo.dstep[idx] = (q2sq != 0.0)
                ? cfg->d0 / q2sq / pow(1.0 + geo.Y / q2sq, 0.13) : 0.0;
        }
    }
    return geo;
}

void geometry_free(Geometry *geo) {
    free(geo->SN); free(geo->Q); free(geo->invQ);
    free(geo->dstep); free(geo->wrap);
    geo->SN = geo->Q = geo->invQ = geo->dstep = NULL;
    geo->wrap = NULL;
}

/* ---- replica lifecycle ----------------------------------------------- */
void replica_alloc(Replica *rep, const Geometry *geo) {
    int L = geo->L;
    rep->L = L;
    rep->h  = calloc((size_t)L * L, sizeof(double complex));
    rep->S  = calloc((size_t)L * L, sizeof(double complex));
    rep->dS = calloc((size_t)L * L, sizeof(double complex));
    rep->g  = calloc((size_t)L * L, sizeof(double));
}

void replica_free(Replica *rep) {
    free(rep->h); free(rep->S); free(rep->dS); free(rep->g);
    rep->h = NULL; rep->S = NULL; rep->dS = NULL; rep->g = NULL;
}

void replica_init(Replica *rep, const Geometry *geo, uint64_t seed, uint64_t stream) {
    int L = geo->L;
    pcg32_seed(&rep->rng, seed, stream);
    /* Start from the harmonic ground-state ensemble h_q = 1/q^2 (real). */
    for (int q1 = 0; q1 < L; q1++)
        for (int q2 = 0; q2 < L; q2++) {
            int i = IDX(q1, q2, L);
            rep->h[i] = (q1 == 0 && q2 == 0) ? (1.0 / (geo->a * geo->a))
                                             : (1.0 / geo->Q[i]);
            rep->g[i] = 0.0;
        }
    rep->nmeas = 0;
    rep->sum_Kx = rep->sum_Ky = rep->sum_KxKx = rep->sum_KxKy = 0.0;
    rep->proposed = rep->accepted = 0;
    calcS(rep, geo);
}

/* ---- full nonlinear array (O(L^4)) ----------------------------------- */
void calcS(Replica *rep, const Geometry *geo) {
    int L = geo->L;
    const double *SN = geo->SN;
    const double complex *h = rep->h;
    const int *wrap = geo->wrap;
    for (int q1 = 0; q1 < L; q1++)
        for (int q2 = 0; q2 < L; q2++) {
            int qi = IDX(q1, q2, L);
            rep->S[qi] = 0.0;
            if (q1 == 0 && q2 == 0) continue;
            double complex acc = 0.0;
            double invQ = geo->invQ[qi];
            for (int k1 = 0; k1 < L; k1++)
                for (int k2 = 0; k2 < L; k2++) {
                    double p = SN[k1] * SN[q2] - SN[k2] * SN[q1];
                    p *= p * invQ;
                    int a1 = wrap[k1 + q1 + L], a2 = wrap[k2 + q2 + L];
                    acc += p * h[IDX(k1, k2, L)] * conj(h[IDX(a1, a2, L)]);
                }
            rep->S[qi] = acc;
        }
}

/* ---- one Metropolis single-mode update ------------------------------- */
static void metropolis_step(Replica *rep, const Geometry *geo, const Config *cfg) {
    int L = geo->L, n = cfg->n;
    const double *SN = geo->SN;
    const double *invQ = geo->invQ;
    const int *wrap = geo->wrap;
    double complex *h = rep->h, *S = rep->S, *dS = rep->dS;

    /* Pick a mode inside the effective zone [-n, n]^2, excluding (0,0). */
    int k1 = (int)(pcg32_next(&rep->rng) % (2u * n + 1u)) - n;
    int k2 = (int)(pcg32_next(&rep->rng) % (2u * n + 1u)) - n;
    if (k1 == 0 && k2 == 0) return;
    if (k1 < 0) k1 += L;
    if (k2 < 0) k2 += L;

    int ki = IDX(k1, k2, L);
    double A = geo->Q[ki];
    double d = geo->dstep[ki];
    double A2 = A * A;

    double complex z = (pcg32_double(&rep->rng) - 0.5) +
                       (pcg32_double(&rep->rng) - 0.5) * I;

    double w = 0.0;
    for (int q1 = 0; q1 < L; q1++)
        for (int q2 = 0; q2 < L; q2++) {
            int qi = IDX(q1, q2, L);
            if (q1 == 0 && q2 == 0) continue;

            double p = SN[k1] * SN[q2] - SN[k2] * SN[q1];  p *= p;
            double kq = SN[wrap[k2 + q2 + L]] * SN[q1] -
                        SN[wrap[k1 + q1 + L]] * SN[q2];      kq *= kq;
            double qk = SN[wrap[L + k1 - q1]] * SN[q2] -
                        SN[wrap[L + k2 - q2]] * SN[q1];      qk *= qk;

            double complex s =
                  p  * conj(h[IDX(wrap[k1 + q1 + L], wrap[k2 + q2 + L], L)]) * z
                + kq * h[IDX(wrap[2 * L - k1 - q1], wrap[2 * L - k2 - q2], L)] * z
                + qk * h[IDX(wrap[k1 - q1 + L], wrap[k2 - q2 + L], L)] * conj(z)
                + p  * conj(h[IDX(wrap[q1 - k1 + L], wrap[q2 - k2 + L], L)]) * conj(z);
            if ((q1 + 2 * k1) % L == 0 && (q2 + 2 * k2) % L == 0)
                s += p * z * z * d;
            if (((q1 - 2 * k1) % L + L) % L == 0 && ((q2 - 2 * k2) % L + L) % L == 0)
                s += p * conj(z) * conj(z) * d;
            s *= d * invQ[qi];
            dS[qi] = s;
            w += creal((2.0 * S[qi] + s) * conj(s));
        }
    w *= -geo->Y / L / L;
    w -= A2 * creal((2.0 * h[ki] + d * z) * conj(z)) * d;

    rep->proposed++;
    if (w > log(pcg32_double(&rep->rng) + 1e-300)) {
        rep->accepted++;
        h[ki] += d * z;
        h[IDX(wrap[L - k1], wrap[L - k2], L)] += d * conj(z);
        for (int i = 0; i < L * L; i++) S[i] += dS[i];
    }
}

void replica_sweep(Replica *rep, const Geometry *geo, const Config *cfg) {
    int l = 2 * cfg->n + 1;
    long steps = (long)l * l;
    for (long s = 0; s < steps; s++) metropolis_step(rep, geo, cfg);
}

/* ---- observables ----------------------------------------------------- */
void replica_measure(Replica *rep, const Geometry *geo) {
    int L = geo->L, N8 = geo->N8;
    const double *SN = geo->SN;
    const double complex *h = rep->h;

    for (int i = 0; i < L * L; i++) rep->g[i] += creal(h[i] * conj(h[i]));
    rep->nmeas++;

    /* Poisson-ratio correlators over the small-q window |q| < q8. */
    double Kx = 0.0, Ky = 0.0;
    for (int q1 = -N8; q1 <= N8; q1++)
        for (int q2 = -N8; q2 <= N8; q2++) {
            if (q1 == 0 && q2 == 0) continue;
            int i1 = geo->wrap[q1 + L], i2 = geo->wrap[q2 + L];
            double hh = creal(h[IDX(i1, i2, L)] * conj(h[IDX(i1, i2, L)]));
            Kx += SN[i1] * SN[i1] * hh;
            Ky += SN[i2] * SN[i2] * hh;
        }
    rep->sum_Kx  += Kx;
    rep->sum_Ky  += Ky;
    rep->sum_KxKx += Kx * Kx;
    rep->sum_KxKy += Kx * Ky;
}

/* ---- reduction across replicas --------------------------------------- */
Result result_reduce(const Replica *reps, int nreps, const Geometry *geo) {
    int L = geo->L;
    Result res;
    res.N = geo->N; res.L = L; res.N8 = geo->N8;
    res.a = geo->a; res.p8 = sqrt(geo->Y * 3.0 / (2.0 * M_PI));
    res.G = calloc((size_t)L * L, sizeof(double));

    long total_meas = 0, proposed = 0, accepted = 0;
    double sKx = 0, sKy = 0, sKxKx = 0, sKxKy = 0;
    for (int r = 0; r < nreps; r++) {
        for (int i = 0; i < L * L; i++) res.G[i] += reps[r].g[i];
        total_meas += reps[r].nmeas;
        proposed  += reps[r].proposed;
        accepted  += reps[r].accepted;
        sKx += reps[r].sum_Kx; sKy += reps[r].sum_Ky;
        sKxKx += reps[r].sum_KxKx; sKxKy += reps[r].sum_KxKy;
    }
    if (total_meas > 0)
        for (int i = 0; i < L * L; i++) res.G[i] /= total_meas;
    res.total_meas = total_meas;
    res.accept_rate = proposed ? (double)accepted / proposed : 0.0;

    /* nu = -(<Kx Ky> - <Kx><Ky>) / (<Kx Kx> - <Kx>^2) */
    res.poisson = 0.0;
    if (geo->N8 > 0 && total_meas > 0) {
        double mKx = sKx / total_meas, mKy = sKy / total_meas;
        double cov = sKxKy / total_meas - mKx * mKy;
        double var = sKxKx / total_meas - mKx * mKx;
        if (var != 0.0) res.poisson = -cov / var;
    }
    return res;
}

void result_free(Result *res) { free(res->G); res->G = NULL; }
