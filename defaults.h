#define nt0 2
#define N0 30
#define n0 30
#define MT0 400
#define M0 1600
#define d0 2.
#define p80 .32    // changes from 0 to pi

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#define pi M_PI

#define rank omp_get_thread_num()
#define NT omp_get_num_threads()

#define PRINT_NOW {time_t rawtime; time(&rawtime); printf("\nnow =\t%s", asctime(localtime(&rawtime)));}
#define SHOW_TIME(i,M) {double tt = omp_get_wtime()-t0, pc = 1.*(i+1)/M;\
        printf("  %.1lf min  \t%.0lf%%\test: %.1f min\tleft: %.1f min\t \r", tt/60, 100.*pc, tt/pc/60, (1-pc)*(tt/pc)/60); fflush(stdout);}


int init(int N, double complex **h, double complex **S, double **g, double **g2, int ***c, double *SN, double **Q);
int dump(int N, double complex **h, double complex **S, double **g, double **g2, int ***c);
int simulate(double p8, int N, int n, double complex **h, double complex **S, double **g, double **g2, int ***c, double *SN, double **Q);

