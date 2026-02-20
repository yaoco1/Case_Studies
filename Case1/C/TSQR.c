#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include "TSQR.h"

void matrix_multiply(double *A, double *B, double *C, int m, int k, int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i*n + j] = 0.0;
            for (int t = 0; t < k; t++)
                C[i*n + j] += A[i*k + t] * B[t*n + j];
        }
    }
}

void TSQR(double *A, int m, int n, double *Q, double *R)
{
    int block_sz = m / 4;
    double *Q_local[4];
    double *R_local[4];
    
    /* local QR */
    for (int b = 0; b < 4; b++) {
        Q_local[b] = malloc(block_sz*n*sizeof(double));
        R_local[b] = malloc(n*n*sizeof(double));

        memcpy(Q_local[b], A + b*block_sz*n, block_sz*n*sizeof(double));

        double *tau = malloc(n*sizeof(double));

        LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, block_sz, n, Q_local[b], n, tau);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                R_local[b][i*n+j] = (j >= i) ? Q_local[b][i*n+j] : 0.0;
            }
        }

        LAPACKE_dorgqr(LAPACK_ROW_MAJOR, block_sz, n, n, Q_local[b], n, tau);

        free(tau);
    }
    
    /* First tree level */
    double *R12 = malloc(2*n*n*sizeof(double));
    double *R34 = malloc(2*n*n*sizeof(double));

    memcpy(R12, R_local[0], n*n*sizeof(double));
    memcpy(R12 + n*n, R_local[1], n*n*sizeof(double));
    memcpy(R34, R_local[2], n*n*sizeof(double));
    memcpy(R34 + n*n, R_local[3], n*n*sizeof(double));

    double tau2[n];
    double tau3[n];

    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, 2*n, n, R12, n, tau2);

    double *Q12 = malloc(2*n*n*sizeof(double));
    memcpy(Q12, R12, 2*n*n*sizeof(double));

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, 2*n, n, n, Q12, n, tau2);

    /* Extracting R12 */
    double *R12_1 = malloc(n*n*sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i)
                R12_1[i*n+j] = R12[i*n+j];
            else
                R12_1[i*n+j] = 0.0;
        }
    }

    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, 2*n, n, R34, n, tau3);

    double *Q34 = malloc(2*n*n*sizeof(double));
    memcpy(Q34, R34, 2*n*n*sizeof(double));

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, 2*n, n, n, Q34, n, tau3);

    double *R34_1 = malloc(n*n*sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i)
                R34_1[i*n+j] = R34[i*n+j];
            else
                R34_1[i*n+j] = 0.0;
        }
    }
    
    /* Second tree level */
    double *R_final = malloc(2*n*n*sizeof(double));
    memcpy(R_final, R12_1, n*n*sizeof(double));
    memcpy(R_final + n*n, R34_1, n*n*sizeof(double));

    double tau4[n];

    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, 2*n, n, R_final, n, tau4);

    double *QF = malloc(2*n*n*sizeof(double));
    memcpy(QF, R_final, 2*n*n*sizeof(double));

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, 2*n, n, n, QF, n, tau4);
    
    /* Extracting final R */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i)
                R[i*n+j] = R_final[i*n+j];
            else
                R[i*n+j] = 0.0;
        }
    }
    
    /* Building final Q */
    for (int b = 0; b < 4; b++) {
        double *Q_temp = malloc(block_sz*n*sizeof(double));
        double *Q_sub1 = malloc(n*n*sizeof(double));
        double *Q_sub2 = malloc(n*n*sizeof(double));

        /* ---- First tree level ---- */
        if (b == 0)
            memcpy(Q_sub1, Q12, n*n*sizeof(double));
        else if (b == 1)
            memcpy(Q_sub1, Q12 + n*n, n*n*sizeof(double));
        else if (b == 2)
            memcpy(Q_sub1, Q34, n*n*sizeof(double));
        else
            memcpy(Q_sub1, Q34 + n*n, n*n*sizeof(double));

        /* ---- Second tree level ---- */
        if (b < 2)
            memcpy(Q_sub2, QF, n*n*sizeof(double));
        else
            memcpy(Q_sub2, QF + n*n, n*n*sizeof(double));

        matrix_multiply(Q_local[b], Q_sub1, Q_temp, block_sz, n, n);
        matrix_multiply(Q_temp, Q_sub2, Q + b*block_sz*n, block_sz, n, n);

        free(Q_temp);
        free(Q_sub1);
        free(Q_sub2);
    }
}
