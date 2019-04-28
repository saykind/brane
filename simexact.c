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
    
    double A = a*a*(k1*k1+k2*k2), d = d0/A/pow(1+p8*p8/a/a/(k1*k1+k2*k2),.2); A = A*A;
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    double complex **SE = (double complex **)malloc(L*sizeof(double complex*));                         // exact
    double complex **dSE= (double complex **)malloc(L*sizeof(double complex*));                         // exact
    for (i = 0; i < L; i++) {
        *(dS+i) = (double complex *)malloc(L*sizeof(double complex));
        *(SE+i)  = (double complex *)malloc(L*sizeof(double complex));                                  // exact
        *(dSE+i) = (double complex *)malloc(L*sizeof(double complex));                                  // exact
    }
    dS[0][0] = 0;

    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    #pragma omp parallel for collapse(2) reduction(+:w)
    for (q1 = -N; q1 < N+1; q1++)
        for (q2 = -N; q2 < N+1; q2++) {
            if (!q1 && !q2) continue;
            SE[(q1+L)%L][(q2+L)%L] = 0;                                                                 // exact
            int p1, p2;                                                                                 // exact
            for (p1 = -N; p1 < N+1; p1++)                                                               // exact
                for (p2 = -N; p2 < N+1; p2++) {                                                         // exact
                    SE[(q1+L)%L][(q2+L)%L] += a*a*(p1*q2-p2*q1)*(p1*q2-p2*q1)/(q1*q1+q2*q2)*h[(p1+L)%L][(p2+L)%L]*conj(h[(p1+q1+L)%L][(p2+q2+L)%L]);
                }                                                                                       // exact
            dSE[(q1+L)%L][(q2+L)%L] = -SE[(q1+L)%L][(q2+L)%L];                                          // exact
            int kq1 = (L-k1-q1)%L < N+1 ? (L-k1-q1)%L : (L-k1-q1)%L-L;
            int kq2 = (L-k2-q2)%L < N+1 ? (L-k2-q2)%L : (L-k2-q2)%L-L;
            int qk1 = (L+k1-q1)%L < N+1 ? (L+k1-q1)%L : (L+k1-q1)%L-L;
            int qk2 = (L+k2-q2)%L < N+1 ? (L+k2-q2)%L : (L+k2-q2)%L-L;
            double complex s;
            s  =  (k1*q2-k2*q1) * (k1*q2-k2*q1) *( conj(h[(k1+q1+L)%L][(k2+q2+L)%L])*z );
            s += (kq1*q2-kq2*q1)*(kq1*q2-kq2*q1)*( h[(L-k1-q1)%L][(L-k2-q2)%L]*z   );
            s += (qk1*q2-qk2*q1)*(qk1*q2-qk2*q1)*( h[(k1-q1+L)%L][(k2-q2+L)%L]*conj(z) );
            s +=  (k1*q2-k2*q1) * (k1*q2-k2*q1) *( conj(h[(q1-k1+L)%L][(q2-k2+L)%L])*conj(z) );
            if (!((q1+2*k1)%L) && !((q2+2*k2)%L))
                s += (k1*q2-k2*q1)*(k1*q2-k2*q1)*z*z*d;
            if (!((q1-2*k1+2*L)%L) && !((q2-2*k2+2*L)%L))
                s += (k1*q2-k2*q1)*(k1*q2-k2*q1)*conj(z)*conj(z)*d;
            s *= a*a*d/(q1*q1+q2*q2);
            dS[(q1+L)%L][(q2+L)%L] = s;
            w += creal( (2*S[(q1+L)%L][(q2+L)%L]+s)*conj(s) );
        }
    w *= -Y/L/L;
    w -= A*creal((2*h[(k1+L)%L][(k2+L)%L] + d*z)*conj(z))*d;
    double E = 0;                                                                                       // exact
    #pragma omp parallel for collapse(2) reduction(+:E)                                                 // exact
    for (q1 = 0; q1 < L; q1++)                                                                          // exact
        for (q2 = 0; q2 < L; q2++)                                                                      // exact
            E += .5*A*creal(h[q1][q2]*conj(h[q1][q2]))+(Y/L/L)*creal(SE[q1][q2]*conj(SE[q1][q2]));      // exact
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
        #pragma omp parallel for collapse(2)                                                            // exact
        for (q1 = -N; q1 < N+1; q1++)                                                                   // exact
            for (q2 = -N; q2 < N+1; q2++) {                                                             // exact
                if (!q1 && !q2) continue;                                                               // exact
                SE[(q1+L)%L][(q2+L)%L] = 0;                                                             // exact
                int p1, p2;                                                                             // exact
                for (p1 = -N; p1 < N+1; p1++)                                                           // exact
                    for (p2 = -N; p2 < N+1; p2++) {                                                     // exact
                        SE[(q1+L)%L][(q2+L)%L] += a*a*(p1*q2-p2*q1)*(p1*q2-p2*q1)/(q1*q1+q2*q2)*h[(p1+L)%L][(p2+L)%L]*conj(h[(p1+q1+L)%L][(p2+q2+L)%L]);
                        dSE[(q1+L)%L][(q2+L)%L] += a*a*(p1*q2-p2*q1)*(p1*q2-p2*q1)/(q1*q1+q2*q2)*h[(p1+L)%L][(p2+L)%L]*conj(h[(p1+q1+L)%L][(p2+q2+L)%L]);
                    }                                                                                   // exact
                if (cabs(dSE[(q1+L)%L][(q2+L)%L]-dS[(q1+L)%L][(q2+L)%L]) > 1e-12)                       // exact
                    printf("(%d,%d)[%d,%d] dSE =\t%.7lf+%.7lfi\tdSE-dS =\t%.15lf\n",\
                    k1, k2, q1, q2, creal(dSE[(q1+L)%L][(q2+L)%L]), cimag(dSE[(q1+L)%L][(q2+L)%L]),\
                    cabs(dSE[(q1+L)%L][(q2+L)%L]-dS[(q1+L)%L][(q2+L)%L]));                              // exact
            }                                                                                           // exact
        #pragma omp parallel for collapse(2) reduction(+:E)                                             // exact
        for (q1 = 0; q1 < L; q1++)                                                                      // exact
            for (q2 = 0; q2 < L; q2++)                                                                  // exact
                E -= .5*A*creal(h[q1][q2]*conj(h[q1][q2]))+(Y/L/L)*creal(SE[q1][q2]*conj(SE[q1][q2]));  // exact
        printf("(%d,%d) E =\t%.7lf\tE-w =\t%.15lf\n", k1, k2, E, E- w);                                 // exact
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
    for (i = 0; i < L; i++) {free(*(dS+i)); free(*(SE+i)); free(*(dSE+i));}                             // exact
    free(dS); free(SE); free(dSE);                                                                      // exact
    return 0;
}
