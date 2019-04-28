#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "defaults.h"

int main(int argc, char *argv[]) {
    int i, j, nt=nt0, N=N0, M=M0, MTH=MT0, L=2*N+1, verbose=0;
    double p8=p80, T0=16e-8, t0; srand(time(NULL));
    for(i = 1; i < argc; i++ ) {
        if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
            printf("Usage:\n");
            printf(" -h\tprint this message\n");
            printf(" -v\tverbose mode [not recommended for N<20]\n");
            printf(" -d\tprint default values\n");
            printf(" nt=%%d\tset number of threads\n");
            printf(" N=%%d\tset linear lattice size L=2*N+1 total size L*L\n");
            printf(" MTH=%%d\tset number of thermalization MC steps (total: MTH*L*L) [recommended value MTH=2000]\n");
            printf(" M=%%d\tset number of MC steps (total: M*L*L)\n");
            printf(" p8=%%f\tset interaction force (0 < p8 < 3.14)\n");
            return 0;
        } 
        if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose") ) {
            verbose = 1;
            continue;
        }
        if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--default") ) {
            printf("Default values:\n");
            printf(" N=%d\t L=%d\t M=%d\t MTH=%d thrds=%d\t p8=%.2lf\n", N, L, M, MTH, nt, p8);
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
    L=2*N+1; omp_set_num_threads(nt);  
    printf("\nN =\t%d\nL =\t%d\nthrds =\t%d\np8 =\t%.2lf\n", N, L, nt, p8);

    /* Memory allocation */
    double complex **h = (double complex **)malloc(L*sizeof(double complex*));
    double complex **S = (double complex **)malloc(L*sizeof(double complex*));
    double **g  = (double **)malloc(L*sizeof(double *));
    double **g2 = (double **)malloc(L*sizeof(double *));
    int ***c = malloc(2*sizeof(int**));
    *(c+0) = (int **)malloc(L*sizeof(int*));
    *(c+1) = (int **)malloc(L*sizeof(int*));
    for (i = 0; i < L; i++) {
        *(h+i) = (double complex *)malloc(L*sizeof(double complex));
        *(S+i) = (double complex *)malloc(L*sizeof(double complex));
        *(g+i) = (double *)malloc(L*sizeof(double));
        *(g2+i) = (double *)malloc(L*sizeof(double));
        *(*(c+0)+i) = (int *)malloc(L*sizeof(int));
        *(*(c+1)+i) = (int *)malloc(L*sizeof(int));
    }

    /* Thermalization */
    if (!init(N,h,S,g,g2,c)) {
        t0 = omp_get_wtime();
        PRINT_NOW	    
        printf("MTH =\t%d*%d\test: %.0f/%d min\n", MTH, L*L, T0*MTH*L*L*L*L/60, nt );
        for (i = 0; i < MTH; i++) {
            for (j = 0; j < L*L; j++)
                simulate(p8, N, h, S, NULL, NULL, NULL);
            if (verbose) SHOW_TIME(i,MTH)
        }
        printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
    }
    /* Monte Carlo average */    
    t0 = omp_get_wtime();
    PRINT_NOW
    printf("M =\t%d*%d\test: %.0f/%d min\n", M, L*L, T0*M*L*L*L*L/60, nt);
    for (i = 0; i < M; i++) {    
        for (j = 0; j < L*L; j++)
            simulate(p8, N, h, S, g, g2, c);
        if (verbose) SHOW_TIME(i,M)
    }
    printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
    PRINT_NOW

    /* Data dump */
    dump(N,h,S,g,g2,c);

    /* Memory de-allocation */ 
    for (i = 0; i < L; i++) 
        {free(*(h+i)); free(*(g+i)); free(*(g2+i)); free(*(S+i)); free(*(*(c+0)+i)); free(*(*(c+1)+i));}
    free(h); free(g); free(g2); free(S); free(*(c+0)); free(*(c+1)); free(c);
    return 0;
}
