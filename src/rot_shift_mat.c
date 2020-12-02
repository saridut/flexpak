#include "rot_shift_mat.h"
#include "utils_math.h"
#include <math.h>

/******************************************************************************/

void rotate_vectors (const int n, const double* v, const double* rotmat,
        double* vrot) {

    double rotmat_trans[9]; /* 3 x 3 array as an 1-D array*/ 

    transpose(3, 3, rotmat, rotmat_trans);
    mm_mult(n, 3, v, 3, rotmat_trans, vrot);
}

/******************************************************************************/

void extract_axis_angle_from_rotmat(const double* rotmat, double axis[3],
       double* angle){

    const int n = 3;
    double tr, s;
    double tmp;
    int k;

    tr = trace(3, rotmat);
    *angle = acos(0.5*(tr-1.0));

    if (*angle > 0) {
        if (*angle < M_PI) {
            axis[0] = rotmat[n*2+1] - rotmat[n*1+2];
            axis[1] = rotmat[n*0+2] - rotmat[n*2+0];
            axis[2] = rotmat[n*1+0] - rotmat[n*0+1];
        }
        else {
            /* Find the location of the largest entry in the diagonal of rotmat */
            k = 0; tmp = rotmat[n*0+0];
            if (rotmat[n*1+1] > rotmat[n*0+0]) {
                tmp = rotmat[n*1+1];
                k = 1;
            }
            if (rotmat[n*2+2] > tmp){
                k = 2;
            }

            switch (k) {
                case 0 :
                    axis[0] = sqrt(rotmat[n*0+0]-rotmat[n*1+1]-rotmat[n*2+2]+1)/2;
                    s = 1.0/(2*axis[0]);
                    axis[1] = s*rotmat[n*0+1];
                    axis[2] = s*rotmat[n*0+2];
                    break;
                case 1 :
                    axis[1] = sqrt(rotmat[n*1+1]-rotmat[n*0+0]-rotmat[n*2+2]+1)/2;
                    s = 1.0/(2*axis[1]);
                    axis[0] = s*rotmat[n*0+1];
                    axis[2] = s*rotmat[n*1+2];
                    break;
                case 2 :
                    axis[2] = sqrt(rotmat[n*2+2]-rotmat[n*0+0]-rotmat[n*1+1]+1)/2;
                    s = 1.0/(2*axis[2]);
                    axis[0] = s*rotmat[n*0+2];
                    axis[1] = s*rotmat[n*1+2];
                    break;
            }
        }
    }
    else {
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
    }

    /* Ensuring that angle is [0, pi) */
    *angle = fmod(*angle, 2*M_PI);

    if (*angle < 0.0) {
        *angle = -*angle;
        axis[0] = -axis[0]; axis[1] = -axis[1]; axis[2] = -axis[2];
    }
    else if (*angle > M_PI) {
        *angle = 2*M_PI - *angle;
        axis[0] = -axis[0]; axis[1] = -axis[1]; axis[2] = -axis[2];
    }
}

/******************************************************************************/

void shift_vectors (const int n, const double* v, const double* shiftmat,
        const bool forward, double* vshift){

    double shiftmat_trans[9]; /* 3 x 3 array as an 1-D array*/ 

    transpose(3,3, shiftmat, shiftmat_trans);
    mm_mult(n, 3, v, 3, shiftmat_trans, vshift);
}

/******************************************************************************/

void shift_tensor2 (const double* A, const double* shiftmat,
        const bool forward, double* Ashift){

    const int n=3;

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            Ashift[n*i+j] = 0.0;
        }
    }

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int p=0; p<n; ++p){
                for (int q=0; q<n; ++q){
                    Ashift[n*i+j] += shiftmat[n*i+p]*shiftmat[n*j+q]*A[n*p+q];
                }
            }
        }
    }

}

/******************************************************************************/

void shift_tensor3 (const double* A, const double* shiftmat,
        const bool forward, double* Ashift){

    const int n = 3;

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int k=0; k<n; ++k){
                Ashift[n*n*i+n*j+k] = 0.0;
            }
        }
    }

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int k=0; k<n; ++k){
                for (int p=0; p<n; ++p){
                    for (int q=0; q<n; ++q){
                        for (int r=0; r<n; ++r){
                            Ashift[n*n*i+n*j+k] 
                                += shiftmat[n*i+p]*shiftmat[n*j+q]*shiftmat[n*k+r]
                                *A[n*n*p+n*q+r];
                        }
                    }
                }
            }
        }
    }

}

/******************************************************************************/

