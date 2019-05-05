#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int simulate(double p8, int N, int n, double complex **h, double complex **S, double **g, double **g2, int ***c, double *SN, double **Q) {
    int L=2*N+1, l=2*n+1, i;
    double a = 2*pi/L, Y = (2*pi)/3*p8*p8, w = 0;
    
    int k1 = rand()%l-n, k2 = rand()%l-n, q1, q2;
    if (k1 < 0) k1 += L; if (k2 < 0) k2 += L;
    if (!k1 && !k2) {simulate(p8, N, n, h, S, g, g2, c, SN, Q); return 0;}
    
    double A = Q[k1][k2], d = d0/A/pow(1+sin(p8)*sin(p8)/A,.18); A *= A;
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    for (i = 0; i < L; i++) 
        *(dS+i) = (double complex *)malloc(L*sizeof(double complex));
    dS[0][0] = 0;

    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    #pragma omp parallel for collapse(2) reduction(+:w)
    for (q1 = 0; q1 < L; q1++)
        for (q2 = 0; q2 < L; q2++) {
            if (!q1 && !q2) continue;
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
    }
    if (c && g && g2) {
        a = creal(h[k1][k2]*conj(h[k1][k2]));
        g[k1][k2] += a;
        g[(L-k1)%L][(L-k2)%L] += a;
        g2[k1][k2] += a*a;
        g2[(L-k1)%L][(L-k2)%L] += a*a;
        c[1][k1][k2]++;
        c[1][(L-k1)%L][(L-k2)%L]++;
    }  
    for (i = 0; i < L; i++) free(*(dS+i));
    free(dS);
    return 0;
}
