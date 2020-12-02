#ifndef ROT_SHIFT_MAT_H
#define ROT_SHIFT_MAT_H

#include <stdbool.h>

void rotate_vectors (const int n, const double* v, const double* rotmat,
        double* vrot);

void extract_axis_angle_from_rotmat(const double* rotmat, double axis[3],
       double* angle);

void shift_vectors (const int n, const double* v, const double* shiftmat,
        const bool forward, double* vshift);

void shift_tensor2 (const double* A, const double* shiftmat,
        const bool forward, double* Ashift);

void shift_tensor3 (const double* A, const double* shiftmat,
        const bool forward, double* Ashift);

#endif /* ROT_SHIFT_MAT_H */
