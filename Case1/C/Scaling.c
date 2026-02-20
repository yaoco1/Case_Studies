#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "TSQR.h"

double wall_time()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec*1e-9;
}

int main()
{
    printf("--- Scaling m ---\n");
    printf("m     time\n");
    
    int n = 10;

    for (int m = 2000; m <= 20000; m += 2000)
    {
        double *A = malloc(m*n*sizeof(double));
        double *Q = malloc(m*n*sizeof(double));
        double *R = malloc(n*n*sizeof(double));

        for (int i = 0; i < m*n; i++)
            A[i] = (double)rand()/RAND_MAX;

        double t1 = wall_time();
        TSQR(A, m, n, Q, R);
        double t2 = wall_time();

        printf("%d  %f\n", m, t2 - t1);

        free(A);
        free(Q);
        free(R);
    }

    printf("\n--- Scaling n ---\n");
    printf("n   time\n");
    
    int m_fixed = 10000;

    for (int n_test = 5; n_test <= 40; n_test += 5)
    {
        double *A = malloc(m_fixed*n_test*sizeof(double));
        double *Q = malloc(m_fixed*n_test*sizeof(double));
        double *R = malloc(n_test*n_test*sizeof(double));

        for (int i = 0; i < m_fixed*n_test; i++)
            A[i] = (double)rand()/RAND_MAX;

        double t1 = wall_time();
        TSQR(A, m_fixed, n_test, Q, R);
        double t2 = wall_time();

        printf("%d   %f\n", n_test, t2 - t1);

        free(A);
        free(Q);
        free(R);
    }

    return 0;
}
