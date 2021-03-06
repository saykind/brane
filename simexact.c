#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int simulate(double Y, int N, int n, double complex **h, double complex **S, double complex **dS, double **g, int ***c, double *SN, double **Q) {
    int L = 2*N+1, l = 2*n+1;    
    int k1 = rand()%l-n, k2 = rand()%l-n, q1, q2;
    if (!k1 && !k2) {simulate(Y, N, n, h, S, dS, g, c, SN, Q); return 0;}
    if (k1 < 0) k1 += L; if (k2 < 0) k2 += L;
    double w = 0, A = Q[k1][k2], d = d0/A/pow(1+Y/A,.13); A *= A;

    double complex **SE = (double complex **)malloc(L*sizeof(double complex*));                         // exact
    double complex **dSE= (double complex **)malloc(L*sizeof(double complex*));                         // exact
    for (q1 = 0; q1 < L; q1++) {                                                                        // exact
        *(SE+q1)  = (double complex *)malloc(L*sizeof(double complex));                                 // exact
        *(dSE+q1) = (double complex *)malloc(L*sizeof(double complex));                                 // exact
    }                                                                                                   // exact

    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    #pragma omp parallel for collapse(2) reduction(+:w)
    for (q1 = 0; q1 < L; q1++)
        for (q2 = 0; q2 < L; q2++) {
            if (!q1 && !q2) continue;
            SE[q1][q2] = 0;                                                                             // exact
            int p1, p2;                                                                                 // exact
            for (p1 = 0; p1 < L; p1++)                                                                  // exact
                for (p2 = 0; p2 < L; p2++) {                                                            // exact
                    double p = SN[p1]*SN[q2]-SN[p2]*SN[q1]; p *= p;                                     // exact
                    p /= Q[q1][q2];                                                                     // exact
                    SE[q1][q2] += p*h[p1][p2]*conj(h[(p1+q1)%L][(p2+q2)%L]);                            // exact
                }                                                                                       // exact
            dSE[q1][q2] = -SE[q1][q2];                                                                  // exact
            double p = SN[k1]*SN[q2]-SN[k2]*SN[q1]; p *= p;
            double kq = SN[(k2+q2)%L]*SN[q1]-SN[(k1+q1)%L]*SN[q2]; kq *= kq;
            double qk = SN[(L+k1-q1)%L]*SN[q2]-SN[(L+k2-q2)%L]*SN[q1]; qk *= qk;
            double complex s = p*conj(h[(k1+q1)%L][(k2+q2)%L])*z;
            s += kq*h[(2*L-k1-q1)%L][(2*L-k2-q2)%L]*z;
            s += qk*h[(k1-q1+L)%L][(k2-q2+L)%L]*conj(z);
            s +=  p*conj(h[(q1-k1+L)%L][(q2-k2+L)%L])*conj(z);
            if (!((q1+2*k1)%L) && !((q2+2*k2)%L)) s += p*z*z*d;
            if (!((q1-2*k1+2*L)%L) && !((q2-2*k2+2*L)%L)) s += p*conj(z)*conj(z)*d;
            s *= d/Q[q1][q2];
            dS[q1][q2] = s;
            w += creal((2*S[q1][q2]+s)*conj(s));
        }
    w *= -Y/L/L;
    w -= A*creal((2*h[k1][k2] + d*z)*conj(z))*d;
    double E = 0;                                                                                       // exact
    #pragma omp parallel for collapse(2) reduction(+:E)                                                 // exact
    for (q1 = 0; q1 < L; q1++)                                                                          // exact
        for (q2 = 0; q2 < L; q2++)                                                                      // exact
            E += .5*A*creal(h[q1][q2]*conj(h[q1][q2]))+(Y/L/L)*creal(SE[q1][q2]*conj(SE[q1][q2]));      // exact
    if (w > log(1.*rand()/RAND_MAX)) {
        h[k1][k2] += d*z;
        h[(L-k1)%L][(L-k2)%L] += d*conj(z);
        #pragma omp parallel for collapse(2)
        for (q1 = 0; q1 < L; q1++)
            for (q2 = 0; q2 < L; q2++)
                S[q1][q2] += dS[q1][q2];
        if (c) {
            c[0][k1][k2]++;
            c[0][(L-k1)%L][(L-k2)%L]++;
        }
        #pragma omp parallel for collapse(2)                                                            // exact
        for (q1 = 0; q1 < L; q1++)                                                                      // exact
            for (q2 = 0; q2 < L; q2++) {                                                                // exact
                if (!q1 && !q2) continue;                                                               // exact
                SE[q1][q2] = 0;                                                                         // exact
                int p1, p2;                                                                             // exact
                for (p1 = 0; p1 < L; p1++)                                                              // exact
                    for (p2 = 0; p2 < L; p2++) {                                                        // exact
                        double p = SN[p1]*SN[q2]-SN[p2]*SN[q1]; p *= p;                                 // exact
                        p /= Q[q1][q2];                                                                 // exact
                        SE[q1][q2] += p*h[p1][p2]*conj(h[(p1+q1)%L][(p2+q2)%L]);                        // exact
                        dSE[q1][q2] += p*h[p1][p2]*conj(h[(p1+q1)%L][(p2+q2)%L]);                       // exact
                    }                                                                                   // exact
                if (cabs(dSE[q1][q2]-dS[q1][q2]) > 5e-11)                                               // exact
                    printf("(%d,%d)[%d,%d]\tdSE =\t%.7lf+%.7lfi\tdSE-dS =\t%.15lf\n", k1, k2, q1, q2,\
                    creal(dSE[q1][q2]), cimag(dSE[q1][q2]), cabs(dSE[q1][q2]-dS[q1][q2]));              // exact
            }                                                                                           // exact
        #pragma omp parallel for collapse(2) reduction(+:E)                                             // exact
        for (q1 = 0; q1 < L; q1++)                                                                      // exact
            for (q2 = 0; q2 < L; q2++)                                                                  // exact
                E -= .5*A*creal(h[q1][q2]*conj(h[q1][q2]))+(Y/L/L)*creal(SE[q1][q2]*conj(SE[q1][q2]));  // exact
        printf("(%d,%d)\tE =\t%.7lf\tE-w =\t%.15lf\n", k1, k2, E, E- w);                                // exact
    }
    if (c && g) {
        double a = creal(h[k1][k2]*conj(h[k1][k2]));
        g[k1][k2] += a;
        g[(L-k1)%L][(L-k2)%L] += a;
        c[1][k1][k2]++;
        c[1][(L-k1)%L][(L-k2)%L]++;
    }
    for (q1 = 0; q1 < L; q1++) {free(*(SE+q1)); free(*(dSE+q1));}                                       // exact
    free(SE); free(dSE);                                                                                // exact
    return 0;
}

