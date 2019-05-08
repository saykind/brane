#define nt0 2
#define N0 40
#define MT0 300
#define M0 700
#define d0 2.1
#define p80 .3    // effectively changes from 0 to pi

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#define pi M_PI

#define rank omp_get_thread_num()
#define NT omp_get_num_threads()

#define PRINT_NOW {time_t rawtime; time(&rawtime); printf("\nnow =\t%s", asctime(localtime(&rawtime)));}
#define SHOW_TIME(i,M,t0) {double tt = omp_get_wtime()-t0, pc = 1.*(i+1)/M;\
        printf("  %.1lf min  \t%.0lf%%\test: %.1f min\tleft: %.1f min\t \r", tt/60, 100.*pc, tt/pc/60, (1-pc)*(tt/pc)/60); fflush(stdout);}


int init(int N, double complex **h, double **g, int ***c, int *C, double *px, double *SN, double **Q);
int calcS(int N, double complex **h, double complex **S, double *SN, double **Q);
double calcPR(int N, int N8, double **g, int ***c, int C, double *px, double *SN);
int dump(int N, int N8, double complex **h, double complex **S, double **g, int ***c, int C, double *px, double *SN);
int simulate(double Y, int N, int n, double complex **h, double complex **S, double complex **dS, double **g, int ***c, double *SN, double **Q);

