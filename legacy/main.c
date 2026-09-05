#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int main(int argc, char *argv[]) {
    int i, j, q1, q2, C=0, nt=nt0, N=N0, M=M0, n=N0, N8=0, MTH=MT0, L=2*N+1, l=2*n+1, verbose=0;
    double a = 2.*pi/L, p8=p80, Y=(2*pi)/3*p8*p8, T0=1e-7, t0; srand(time(NULL));
    for(i = 1; i < argc; i++ ) {
        if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
            printf("Usage:\n");
            printf(" -h\tprint this message\n");
            printf(" -v\tverbose mode (not recommended for N<30)\n");
            printf(" -d\tprint default values\n");
            printf(" nt=%%d\tset number of threads\n");
            printf(" N=%%d\tset linear lattice size L=2*N+1 (total: L*L)\n");
            printf(" n=%%d\tset linear effective lattice size l=2*n+1 (default: n=N)\n");
            printf(" MTH=%%d\tset number of thermalization MC steps (total: MTH*l*l)\n");
            printf(" M=%%d\tset number of MC steps (total: M*l*l)\n");
            printf(" p8=%%f\tset interaction force (0 < p8 < 3.14)\n");
            return 0;
        } 
        if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose") ) {
            verbose = 1;
            continue;
        }
        if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--default") ) {
            printf("Default values:\n");
            printf(" N=%d\t n=%d\t L=%d\t l=%d\t M=%d\t MTH=%d thrds=%d\t p8=%.2lf\n", N, n, L, l, M, MTH, nt, p8);
            return 0;
        }  
        if (argv[i][0] == '-') {
            printf("Unknown option '%s'\n", argv[i]);
            return 0;
        }           
        sscanf(argv[i], "nt=%d", &nt);
        sscanf(argv[i], "N=%d", &N);
        sscanf(argv[i], "M=%d", &M);
        sscanf(argv[i], "MTH=%d", &MTH);
        sscanf(argv[i], "p8=%lf", &p8);
    }
    omp_set_num_threads(nt);
    n=N; L=2*N+1; l=2*n+1; a = 2.*pi/L; Y=(2*pi)/3*p8*p8;
    for(i = 1; i < argc; i++ ) {sscanf(argv[i], "n=%d", &n);}
    l=2*n+1; N8=(int)(p8/pi*N);
    printf("\nN =\t%d\nn =\t%d\nN8 =\t%d\nL =\t%d\nl =\t%d\nthrds =\t%d\np8 =\t%.2lf\n", N, n, N8, L, l, nt, p8);

    /* Memory allocation */
    double complex **h = (double complex **)malloc(L*sizeof(double complex*));
    double complex **S = (double complex **)malloc(L*sizeof(double complex*));
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    double **g  = (double **)malloc(L*sizeof(double *));
    int ***c = malloc(2*sizeof(int**));
    *(c+0) = (int **)malloc(L*sizeof(int*));
    *(c+1) = (int **)malloc(L*sizeof(int*));
    double *SN  = (double *)malloc(L*sizeof(double));
    double **Q  = (double **)malloc(L*sizeof(double *));
    double *px  = (double *)malloc(2*sizeof(double));
    for (i = 0; i < L; i++) {
        *(h+i) = (double complex *)malloc(L*sizeof(double complex));
        *(S+i) = (double complex *)malloc(L*sizeof(double complex));
        *(dS+i) = (double complex *)malloc(L*sizeof(double complex));
        *(g+i) = (double *)malloc(L*sizeof(double));
        *(*(c+0)+i) = (int *)malloc(L*sizeof(int));
        *(*(c+1)+i) = (int *)malloc(L*sizeof(int));
        *(Q+i) = (double *)malloc(L*sizeof(double));
    }
    dS[0][0] = 0;
    for (q1 = -N; q1 < N+1; q1++) {
        SN[(q1+L)%L] = sin(a*q1);
        for (q2 = -N; q2 < N+1; q2++)
            Q[(q1+L)%L][(q2+L)%L] = 4*(sin(a*q1/2)*sin(a*q1/2)+sin(a*q2/2)*sin(a*q2/2));
    }

    /* Thermalization */
    if (!init(N, h, g, c, &C, px, SN, Q)) {
        calcS(N, h, S, SN, Q);
        t0 = omp_get_wtime();
        PRINT_NOW	    
        printf("MTH =\t%d*%d\test: %.0f min\n", MTH, l*l, T0*MTH*l*l*L*L/60/nt );
        for (i = 0; i < MTH; i++) {
            for (j = 0; j < l*l; j++)
                simulate(Y, N, n, h, S, dS, NULL, NULL, SN, Q);
            if (verbose) SHOW_TIME(i,MTH,t0)
        }
        printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
    }
    /* Monte Carlo average */ 
    if (M) calcS(N, h, S, SN, Q);   
    t0 = omp_get_wtime();
    PRINT_NOW
    printf("M =\t%d*%d\test: %.0f min\n", M, l*l, T0*M*l*l*L*L/60/nt);
    for (i = 0; i < M; i++) {    
        for (j = 0; j < l*l; j++)
            simulate(Y, N, n, h, S, dS, g, c, SN, Q);
        if (verbose) SHOW_TIME(i,M,t0)
        /* Poisson ratio calculation */
        double Kx = 0, Ky = 0, hh = 0;
        #pragma omp parallel for collapse(2) private(hh) reduction(+:Kx,Ky)
        for (q1 = -N8; q1 < N8+1; q1++) 
            for (q2 = -N8; q2 < N8+1; q2++) {
                hh = creal(creal(h[(q1+L)%L][(q2+L)%L]*conj(h[(q1+L)%L][(q2+L)%L])));
                Kx += SN[(q1+L)%L]*SN[(q1+L)%L]*hh;
                Ky += SN[(q2+L)%L]*SN[(q2+L)%L]*hh;
            }
        px[0] += Kx*Kx; px[1] += Kx*Ky;
    }
    C += M;
    printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
    PRINT_NOW

    /* Data dump */
    dump(N, N8, h, S, g, c, C, px, SN);

    /* Memory de-allocation */ 
    for (i = 0; i < L; i++) 
        {free(*(h+i)); free(*(S+i)); free(*(dS+i)); free(*(g+i)); free(*(*(c+0)+i)); free(*(*(c+1)+i)); free(*(Q+i));}
    free(h); free(S); free(dS); free(g); free(*(c+0)); free(*(c+1)); free(c); free(SN); free(Q); free(px);
    return 0;
}
