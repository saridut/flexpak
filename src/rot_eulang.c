#include "rot_eulang.h"
#include "utils_math.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>


/******************************************************************************/

static void rotmat_XYZ(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = sin_phi*sin_theta*cos_psi - cos_phi*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta*cos_psi + sin_phi*sin_psi;
    rotmat[n*1+0] = cos_theta*sin_psi;
    rotmat[n*1+1] = sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
    rotmat[n*1+2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_psi;
    rotmat[n*2+0] = -sin_theta;
    rotmat[n*2+1] = sin_phi*cos_theta;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void rotmat_XZY(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = sin_phi*sin_theta - cos_phi*cos_theta*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta + sin_phi*cos_theta*sin_psi;
    rotmat[n*1+0] = sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = -sin_phi*cos_psi;
    rotmat[n*2+0] = -sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_theta + cos_phi*sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta - sin_phi*sin_theta*sin_psi;
}

/******************************************************************************/

static void rotmat_YXZ(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi - sin_phi*sin_theta*sin_psi;
    rotmat[n*0+1] = -cos_phi*sin_psi;
    rotmat[n*0+2] = sin_theta*cos_psi + sin_phi*cos_theta*sin_psi;
    rotmat[n*1+0] = sin_phi*sin_theta*cos_psi + cos_theta*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = sin_theta*sin_psi - sin_phi*cos_theta*cos_psi;
    rotmat[n*2+0] = -cos_phi*sin_theta;
    rotmat[n*2+1] = sin_phi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void rotmat_YZX(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = -sin_psi;
    rotmat[n*0+2] = sin_theta*cos_psi;
    rotmat[n*1+0] = sin_phi*sin_theta + cos_phi*cos_theta*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_theta;
    rotmat[n*2+0] = sin_phi*cos_theta*sin_psi - cos_phi*sin_theta;
    rotmat[n*2+1] = sin_phi*cos_psi;
    rotmat[n*2+2] = sin_phi*sin_theta*sin_psi + cos_phi*cos_theta;
}

/******************************************************************************/

static void rotmat_ZXY(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi + sin_phi*sin_theta*sin_psi;
    rotmat[n*0+1] = sin_phi*sin_theta*cos_psi - cos_theta*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta;
    rotmat[n*1+0] = cos_phi*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = -sin_phi;
    rotmat[n*2+0] = sin_phi*cos_theta*sin_psi - sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_theta*cos_psi + sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void rotmat_ZYX(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = -cos_theta*sin_psi;
    rotmat[n*0+2] = sin_theta;
    rotmat[n*1+0] = sin_phi*sin_theta*cos_psi + cos_phi*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi - sin_phi*sin_theta*sin_psi;
    rotmat[n*1+2] = -sin_phi*cos_theta;
    rotmat[n*2+0] = sin_phi*sin_psi - cos_phi*sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_psi + cos_phi*sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_XYZ(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*2+0] < 1.0){
        if (rotmat[n*2+0] > -1.0){
            theta = asin(-rotmat[n*2+0]);
            psi = atan2(rotmat[n*1+0], rotmat[n*0+0]);
            phi = atan2(rotmat[n*2+1], rotmat[n*2+2]);
        }
        else {
            /* Not unique: phi - psi = atan2(-rotmat[1][2], rotmat[1][1]) */
            theta = M_PI_2;
            psi = -atan2(-rotmat[n*1+2], rotmat[n*1+1]);
            phi = 0.0;
        }
    }
    else {
        /* Not unique: phi + psi = atan2(-rotmat[1][2], rotmat[1][1]) */
        phi = 0.0;
        theta = -M_PI_2;
        psi = atan2(-rotmat[n*1+2], rotmat[n*1+1]);
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_XZY(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*1+0] < 1.0){
        if (rotmat[n*1+0] > -1.0){
            phi = atan2(-rotmat[n*1+2], rotmat[n*1+1]);
            theta = atan2(-rotmat[n*2+0], rotmat[n*0+0]);
            psi = asin(rotmat[n*1+0]);
        }
        else {
            /* Not unique: phi - theta = atan2(rotmat[2][1], rotmat[2][2]) */
            phi = 0.0;
            theta = -atan2(rotmat[n*2+1], rotmat[n*2+2]);
            psi = -M_PI_2;
        }
    }
    else {
        /* Not unique: phi + theta = atan2(rotmat[2][1], rotmat[2][2]) */
        phi = 0.0;
        theta = atan2(rotmat[n*2+1], rotmat[n*2+2]);
        psi = M_PI_2;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_YXZ(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*2+1] < 1.0) {
        if (rotmat[n*2+1] > -1.0) {
            phi = asin(rotmat[n*2+1]);
            theta = atan2(-rotmat[n*2+0], rotmat[n*2+2]);
            psi = atan2(-rotmat[n*0+1], rotmat[n*1+1]);
        }
        else {
            /* Not unique: theta - psi = atan2(rotmat[0][2], rotmat[0][0]) */
            phi = -M_PI_2;
            theta = 0.0;
            psi = -atan2(rotmat[n*0+2], rotmat[n*0+0]);
        }
    }
    else {
        /* Not unique: theta + psi = atan2(rotmat[0][2], rotmat[0][0]) */
        phi = M_PI_2;
        theta = 0.0;
        psi = atan2(rotmat[n*0+2], rotmat[n*0+0]);
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_YZX(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*0+1] < 1.0) {
        if (rotmat[n*0+1] > -1.0) {
            phi = atan2(rotmat[n*2+1], rotmat[n*1+1]);
            theta = atan2(rotmat[n*0+2], rotmat[n*0+0]);
            psi = asin(-rotmat[n*0+1]);
        }
        else {
            /* Not unique: theta - phi = atan2(-rotmat[2][0], rotmat[2][2]) */
            phi = -atan2(-rotmat[n*2+0], rotmat[n*2+2]);
            theta = 0.0;
            psi = M_PI_2;
        }
    }
    else {
        /* Not unique: theta + phi = atan2(-rotmat[2][0], rotmat[2][2]) */
        phi = atan2(-rotmat[n*2+0], rotmat[n*2+2]);
        theta = 0.0;
        psi = -M_PI_2;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_ZXY(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*1+2] < 1.0) {
        if (rotmat[n*1+2] > -1.0) {
            phi = asin(-rotmat[n*1+2]);
            theta = atan2(rotmat[n*0+2], rotmat[n*2+2]);
            psi = atan2(rotmat[n*1+0], rotmat[n*1+1]);
        }
        else {
            /* Not unique: psi - theta = atan2(-rotmat[0][1], rotmat[0][0]) */
            phi = M_PI_2;
            theta = -atan2(-rotmat[n*0+1], rotmat[n*0+0]);
            psi = 0.0;
        }
    }
    else {
        /* Not unique: psi + theta = atan2(-rotmat[0][1], rotmat[0][0]) */
        phi = -M_PI_2;
        theta = atan2(-rotmat[n*0+1], rotmat[n*0+0]);
        psi = 0.0;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_ZYX(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*0+2] < 1.0) {
        if (rotmat[n*0+2] > -1.0) {
            phi = atan2(-rotmat[n*1+2], rotmat[n*2+2]);
            theta = asin(rotmat[n*0+2]);
            psi = atan2(-rotmat[n*0+1], rotmat[n*0+0]);
        }
        else {
            /* Not unique: psi - phi = atan2(rotmat[1][0], rotmat[1][1]) */
            phi = -atan2(rotmat[n*1+0], rotmat[n*1+1]);
            theta = -M_PI_2;
            psi = 0.0;
        }
    }
    else {
        /* Not unique: psi + phi = atan2(rotmat[1][0], rotmat[1][1]) */
        phi = atan2(rotmat[n*1+0], rotmat[n*1+1]);
        theta = M_PI_2;
        psi = 0.0;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

void rotmat_euler(const double* euler, enum EulangSeq seq, const bool world,
        double* rotmat){

    double euler_[3];
    double rotmat_[9]; /* 3 x 3 matrix as an 1-D array */

    if (!world){
        for (int i=0; i<3; ++i){
            euler_[i] = -euler[i];
        }
    }
    else {
        for (int i=0; i<3; ++i){
            euler_[i] = euler[i];
        }
    }

    switch (seq){
        case ES_XYZ :
            rotmat_XYZ(euler_, rotmat_);
            break;
        case ES_XZY :
            rotmat_XZY(euler_, rotmat_);
            break;
        case ES_YXZ :
            rotmat_YXZ(euler_, rotmat_);
            break;
        case ES_YZX :
            rotmat_YZX(euler_, rotmat_);
            break;
        case ES_ZXY :
            rotmat_ZXY(euler_, rotmat_);
            break;
        case ES_ZYX :
            rotmat_ZYX(euler_, rotmat_);
            break;
    }

    if (!world) {
        transpose(3, 3, rotmat_, rotmat);
    }
    else {
        memcpy(rotmat, rotmat_, 9*sizeof(double));
    }
}

/******************************************************************************/

void factor_rotmat(const double* rotmat, enum EulangSeq seq,
        const bool world, double* euler){

    double rotmat_[9]; /* 3 x 3 matrix as an 1-D array */

    if (!world) {
        transpose(3, 3, rotmat, rotmat_);
    }
    else {
        memcpy(rotmat_, rotmat, 9*sizeof(double));
    }

    switch (seq){
        case ES_XYZ :
            factor_rotmat_XYZ(rotmat_, euler);
            break;
        case ES_XZY :
            factor_rotmat_XZY(rotmat_, euler);
            break;
        case ES_YXZ :
            factor_rotmat_YXZ(rotmat_, euler);
            break;
        case ES_YZX :
            factor_rotmat_YZX(rotmat_, euler);
            break;
        case ES_ZXY :
            factor_rotmat_ZXY(rotmat_, euler);
            break;
        case ES_ZYX :
            factor_rotmat_ZYX(rotmat_, euler);
            break;
    }

    if (!world){
        for (int i=0; i<3; ++i){
            euler[i] = -euler[i];
        }
    }
}

/******************************************************************************/
