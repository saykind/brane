#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int simulate(double p8, int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {
    int L=2*N+1, i;
    double a = 2.*pi/L, Y = (2*pi)/3*p8*p8, w = 0;
    
    int k1 = rand()%L-N, k2 = rand()%L-N, q1, q2;
    if (!k1 && !k2) {simulate(p8, N, h, S, g, g2, c); return 0;}
    
    double A = 4*(sin(a*k1/2)*sin(a*k1/2)+sin(a*k2/2)*sin(a*k2/2)), d = d0/A/pow(1+sin(p8)*sin(p8)/A,.18); A = A*A;
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    for (i = 0; i < L; i++) 
        *(dS+i) = (double complex *)malloc(L*sizeof(double complex));
    dS[0][0] = 0;

    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    #pragma omp parallel for collapse(2) reduction(+:w)
    for (q1 = -N; q1 < N+1; q1++)
        for (q2 = -N; q2 < N+1; q2++) {
            if (!q1 && !q2) continue;
            int kq1 = (L-k1-q1)%L < N+1 ? (L-k1-q1)%L : (L-k1-q1)%L-L;
            int kq2 = (L-k2-q2)%L < N+1 ? (L-k2-q2)%L : (L-k2-q2)%L-L;
            int qk1 = (L+k1-q1)%L < N+1 ? (L+k1-q1)%L : (L+k1-q1)%L-L;
            int qk2 = (L+k2-q2)%L < N+1 ? (L+k2-q2)%L : (L+k2-q2)%L-L;
            double p = (sin(a*k1)*sin(a*q2)-sin(a*k2)*sin(a*q1))*(sin(a*k1)*sin(a*q2)-sin(a*k2)*sin(a*q1));
            double complex s;
            s  =  p*( conj(h[(k1+q1+L)%L][(k2+q2+L)%L])*z );
            s += (sin(a*kq1)*sin(a*q2)-sin(a*kq2)*sin(a*q1))*(sin(a*kq1)*sin(a*q2)-sin(a*kq2)*sin(a*q1))*( h[(L-k1-q1)%L][(L-k2-q2)%L]*z   );
            s += (sin(a*qk1)*sin(a*q2)-sin(a*qk2)*sin(a*q1))*(sin(a*qk1)*sin(a*q2)-sin(a*qk2)*sin(a*q1))*( h[(k1-q1+L)%L][(k2-q2+L)%L]*conj(z) );
            s +=  p*( conj(h[(q1-k1+L)%L][(q2-k2+L)%L])*conj(z) );
            if (!((q1+2*k1)%L) && !((q2+2*k2)%L))
                s += p*z*z*d;
            if (!((q1-2*k1+2*L)%L) && !((q2-2*k2+2*L)%L))
                s += p*conj(z)*conj(z)*d;
            s *= d/(sin(a*q1)*sin(a*q1)+sin(a*q2)*sin(a*q2));
            dS[(q1+L)%L][(q2+L)%L] = s;
            w += creal( (2*S[(q1+L)%L][(q2+L)%L]+s)*conj(s) );
        }
    w *= -Y/L/L;
    w -= A*creal((2*h[(k1+L)%L][(k2+L)%L] + d*z)*conj(z))*d;
    if (w > log(1.*rand()/RAND_MAX)) {
        h[(k1+L)%L][(k2+L)%L] += d*z;
        h[(L-k1)%L][(L-k2)%L] += d*conj(z);
        #pragma omp parallel for collapse(2)
        for (q1 = 0; q1 < L; q1++)
            for (q2 = 0; q2 < L; q2++)
                S[q1][q2] += dS[q1][q2];
        if (c) {
            c[0][(k1+L)%L][(k2+L)%L]++;
            c[0][(L-k1)%L][(L-k2)%L]++;
        }
    }
    if (c && g && g2) {
        a = creal(h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+L)%L][(k2+L)%L]));
        g[(k1+L)%L][(k2+L)%L] += a;
        g[(L-k1)%L][(L-k2)%L] += a;
        g2[(k1+L)%L][(k2+L)%L] += a*a;
        g2[(L-k1)%L][(L-k2)%L] += a*a;
        c[1][(k1+L)%L][(k2+L)%L]++;
        c[1][(L-k1)%L][(L-k2)%L]++;
    }  
    for (i = 0; i < L; i++) free(*(dS+i));
    free(dS);
    return 0;
}
