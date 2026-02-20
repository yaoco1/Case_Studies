#ifndef TSQR_H
#define TSQR_H

void matrix_multiply(double *A, double *B, double *C, int m, int k, int n);
void TSQR(double *A, int m, int n, double *Q, double *R);

#endif
