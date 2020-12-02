#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

/*!
 *  Returns a pointer to the Gauss-Legendre nodes for a given
 *  value of n, where 2 <= n <= 10
 */
double* gleg_get_nodes(const int n);


/*!
 *  Returns a pointer to the Gauss-Legendre nodes for a given
 *  value of n, where 2 <= n <= 10
 */
double* gleg_get_weights(const int n);


#endif //GAUSS_LEGENDRE_H
