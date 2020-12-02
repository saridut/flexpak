#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <complex.h>

double dot(const int n, const double* a, const double* b);

void cross(const double a[3], const double b[3], double c[3]);

void outer(const int na, const int nb, const double* a,
        const double* b, double* c);

void cross_mat(const double a[3], double* ax);

double norm(const int n, const double* a);

/* Calculates a.(b x c) */
double scalar_triple_product(const double a[3], const double b[3],
        const double c[3]);

/* Calculates a x (b x c) */
void vector_triple_product(const double a[3], const double b[3],
        const double c[3], double d[3]);

/* Calculates determinant of n x n matrix A. Currently handles n = 2, 3, & 4.*/
double det(const int n, const double* A);

/* Kronecker delta function, 0 <= i,j <= 2 */
double kron_delta(const int i, const int j);

/* Returns sum of the elements along the diagonal of an n x n square matrix */
double trace(const int n, const double* A);

/* Check if two floating point numbers of type double are close within */
/* tolerance.                                                          */
bool isclose(const double a, const double b,
        const double rel_tol, const double abs_tol);

/* Check if two arrays of type double are elementwise close within */
/* tolerance.                                                      */
bool allclose(const int n, const double* a, const double* b,
        const double rel_tol, const double abs_tol);

/* Normalizes a vector in-place */
void unitize(const int n, double* a);

/* Generates an identity matrix of size n x n */
void eye(const int n, double* a);

/* Performs a matrix vector multiplication: c_j = A_ij b_j*/
void mv_mult(const int nrows, const int ncols, const double* A,
        const double* b, double* c);

/* Performs a matrix matrix multiplication: C_ij = A_ik B_kj */
void mm_mult(const int nrows_A, const int ncols_A, const double* A,
        const int ncols_B, const double* B, double* C);

/* Transposes a matrix: B = transpose(A) */
void transpose(const int nrows, const int ncols, const double* A, double* B);

/* Adds a square matrix and its transpose in place: Aij = Aij + Aji */
void add_transpose(const int n, double* A);

/* Multiplies a matrix with its transpose: B = A * A^T */
void mult_transpose(const int nrows, const int ncols, const double* A,
        double* B);

/* Calculates the mean of a matrix along dimension dim */
void mat_mean(const int nrows, const int ncols, const double* A,
        const int dim, double* out);

/* Calculates the mean of a vector */
double vec_mean(const int n, const double* v);

/*
 *  @brief Zeros out a memory block of length n
 *  @todo Check whether setting all bytes to zero indicates 0.0 in floating
 *  point.  If not, put in a loop.
 */
void zero_out(const int n, double* v);
void izero_out(const int n, int* v);
void zzero_out(const int n, double complex* v);

#ifdef __cplusplus
}
#endif

#endif // UTILS_MATH_H
