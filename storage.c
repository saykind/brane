#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int init(int N, double complex **h, double complex **S, double **g, double **g2, int ***c, double *SN, double **Q) {
    int L=2*N+1, q1, q2;
    char name[20];
    sprintf(name,"data/N=%d.dat", N);
    FILE *file = fopen(name, "r");
    if (file) {
        double re, im;
        for (q1 = 0; q1 < L; q1++)
            for (q2 = 0; q2 < L; q2++) {
                fscanf(file,"%lf\t%lf\n", &re, &im );
                h[q1][q2] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &re, &im );
                S[q1][q2] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &g[q1][q2], &g2[q1][q2] );
                fscanf(file,"%d\t%d\n", &c[0][q1][q2], &c[1][q1][q2] );
            }        
        fclose(file);
        return 1;
    }
    double a = 2.*pi/L;
    h[0][0] = 1/a/a;
    #pragma omp parallel for private(q1, q2) collapse(2)
    for (q1 = 0; q1 < L; q1++)
         for (q2 = 0; q2 < L; q2++) {
            c[0][q1][q2] = 0;
            c[1][q1][q2] = 1;
            S[q1][q2] = 0;
            if (!q1 && !q2) continue;
            h[q1][q2] = 1./Q[q1][q2];
            g[q1][q2] = h[q1][q2]*h[q1][q2];
            g2[q1][q2] = g[q1][q2]*g[q1][q2];
        }  
    #pragma omp parallel for private(q1, q2) collapse(2)
    for (q1 = 0; q1 < L; q1++)
        for (q2 = 0; q2 < L; q2++) {
            if (!q1 && !q2) continue;
            int k1, k2;
            for (k1 = 0; k1 < L; k1++)
                for (k2 = 0; k2 < L; k2++)   {
                    double p = SN[k1]*SN[q2]-SN[k2]*SN[q1]; p *= p;
                    p /= Q[q1][q2];
                    S[q1][q2] += p*h[k1][k2]*conj(h[(k1+q1)%L][(k2+q2)%L]);
                }
        }
    return 0;
}

int dump(int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {    
    int L=2*N+1, i, q1, q2;
    double a = 2*pi/L;
    char name[20];

    mkdir("data", 0777);
    sprintf(name,"data/N=%d.dat", N); 
    FILE *file = fopen(name, "w"); 
    if (!file) {printf("Cannot save data\n"); return -1;}
    for (q1 = 0; q1 < L; q1++)
        for (q2 = 0; q2 < L; q2++) {
            fprintf(file,"%.15lf\t%.15lf\n", creal(h[q1][q2]), cimag(h[q1][q2]) );
            fprintf(file,"%.15lf\t%.15lf\n", creal(S[q1][q2]), cimag(S[q1][q2]) );
            fprintf(file,"%lf\t%lf\n", g[q1][q2], g2[q1][q2] );
            fprintf(file,"%d\t%d\n", c[0][q1][q2], c[1][q1][q2] );
        }
    fclose(file);

    mkdir("eta", 0777);
    sprintf(name, "eta/N=%d",N);
    file = fopen(name, "w");
    sprintf(name, "eta/fit");
    FILE *fite = fopen(name, "a"); 
    if (!file) {printf("Cannot save data\n"); return -1;}
    for (i = 1; i < N+1; i++) {
        fprintf(file,"%lf\t%lf\t%lf\n", i*a,\
         c[1][0][i]/g[0][i], sqrt(fabs(g2[0][i]-g[0][i]*g[0][i]/c[1][0][i]))/c[1][0][i]/g[0][i]/g[0][i]);
        fprintf(fite,"%lf\t%lf\n", i*a, c[1][0][i]/g[0][i]);
        fprintf(file,"%lf\t%lf\t%lf\n", i*a,\
         c[1][i][0]/g[i][0], sqrt(fabs(g2[i][0]-g[i][0]*g[i][0]/c[1][i][0]))/c[1][i][0]/g[i][0]/g[i][0]);
        fprintf(fite,"%lf\t%lf\n", i*a, c[1][i][0]/g[i][0]);
        fprintf(file,"%lf\t%lf\t%lf\n", i*a*sqrt(2),\
         c[1][i][i]/g[i][i], sqrt(fabs(g2[i][i]-g[i][i]*g[i][i]/c[1][i][i]))/c[1][i][i]/g[i][i]/g[i][i]);
        fprintf(fite,"%lf\t%lf\n", i*a*sqrt(2), c[1][i][i]/g[i][i]);
    }
    fclose(file);

/*
    mkdir("3d", 0777);
    sprintf(name, "3d/N=%d",N);
    file = fopen(name, "w"); 
    if (!file) {printf("Cannot save data\n"); return -1;}
    for (q1 = -N; q1 < N+1; q1++)
        for (q2 = -N; q2 < N+1; q2++) {
            if (!q1 && !q2) continue;
            fprintf(file,"%lf\t%lf\t%lf\n", q1*a, q2*a, c[1][(q1+L)%L][(q2+L)%L]/g[(q1+L)%L][(q2+L)%L] );
        }
    fclose(file);
*/
    return 0;
}
