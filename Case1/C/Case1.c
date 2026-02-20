#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TSQR.h"

double norm(double *A, int m, int n)
{
    double sum = 0.0;
    for (int i = 0; i < m*n; i++)
        sum += A[i]*A[i];
    return sqrt(sum);
}

int main()
{
    int m = 1000;
    int n = 10;

    double *A = malloc(m*n*sizeof(double));
    double *Q = malloc(m*n*sizeof(double));
    double *R = malloc(n*n*sizeof(double));
    double *QR = malloc(m*n*sizeof(double));

    for (int i = 0; i < m*n; i++)
        A[i] = (double)rand()/RAND_MAX;

    TSQR(A, m, n, Q, R);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            QR[i*n+j] = 0.0;
            for (int k = 0; k < n; k++)
                QR[i*n+j] += Q[i*n+k]*R[k*n+j];
        }
    }
    
    double error = 0.0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double diff = A[i*n + j] - QR[i*n + j];
            error += diff * diff;
        }
    }

    error = sqrt(error);

    printf("||A - QR|| = %e\n", error);

    free(A);
    free(Q);
    free(R);
    free(QR);

    return 0;
}
