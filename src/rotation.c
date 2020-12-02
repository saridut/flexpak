#include "rotation.h"
#include "utils_math.h"
#include "ran_num.h"
#include <math.h>
#include <string.h>

/* Euler angle sequence: XYZ (world). First rotation about X, second rotation   */
/* about Y, and the third rotation about Z axis of the world(i.e. fixed) frame. */
/* This is the same as the sequence used in Blender.                            */
/* In contrast, the XYZ sequence is understood in the Aerospace community as:   */
/* First rotation about Z-axis, second rotation about Y-axis, and the third     */
/* rotation about X-axis of the body frame.                                     */

/******************************************************************************/
                            /* AXIS  ANGLE */
/******************************************************************************/

/* Returns angle in [0, pi) */
void fix_axis_angle(double axis[3], double* angle, const bool normalize){

    if (normalize) {
        double nrm = norm(3, axis);
        if (!isclose(nrm, 1.0, 1e-14, 1e-14)) {
            for (int i=0; i<3; ++i) {
                axis[i] /= nrm;
            }
        }
    }

    *angle = fmod(*angle, 2*M_PI);

    if (*angle < 0.0) {
        *angle = -*angle;
        for (int i=0; i<3; ++i) {
            axis[i] = -axis[i];
        }
    }
    else if (*angle > M_PI) {
        *angle = 2*M_PI - *angle;
        for (int i=0; i<3; ++i) {
            axis[i] = -axis[i];
        }
    }
}

/******************************************************************************/

/* Generates a random pair of axis-angle. The axis is a random vector from */
/* the surface of a unit sphere. Algorithm from Allen & Tildesley p. 349.  */
void get_rand_axis_angle(double axis[3], double* angle){

    double zeta1, zeta2, zetasq, rt;
    unsigned int seed = get_random_seed();
    RandomStream stream = init_stream(seed);

    /* Generate angle: A uniform random number from [0.0, 2*pi) */
    duni(stream, 0.0, 2*M_PI, 1, angle);

    while (true) {
        /* Generate two uniform random numbers from [-1, 1) */
        duni(stream, -1.0, 1.0, 1, &zeta1);
        duni(stream, -1.0, 1.0, 1, &zeta2);

        zeta1 = 2*zeta1 - 1.0;
        zeta2 = 2*zeta2 - 1.0;
        zetasq = zeta1*zeta1 + zeta2*zeta2;
        if (zetasq <= 1.0) break;
    }

    rt = sqrt(1.0-zetasq);

    axis[0] = 2*zeta1*rt ;
    axis[1] = 2*zeta2*rt ;
    axis[2] = 1 - 2*zetasq ;

    delete_stream(&stream);
}

/******************************************************************************/

void get_rotmat_axis_angle(const double axis[3], const double angle,
        double* rotmat) {

    const int n = 3;
    double sine, cosine, icos;

    sine = sin(angle);
    cosine = cos(angle);
    icos = 1.0 - cosine;

    rotmat[n*0+0] = axis[0]*axis[0]*icos + cosine;
    rotmat[n*0+1] = axis[0]*axis[1]*icos - axis[2]*sine;
    rotmat[n*0+2] = axis[0]*axis[2]*icos + axis[1]*sine;
    rotmat[n*1+0] = axis[0]*axis[1]*icos + axis[2]*sine;
    rotmat[n*1+1] = axis[1]*axis[1]*icos + cosine;
    rotmat[n*1+2] = axis[1]*axis[2]*icos - axis[0]*sine;
    rotmat[n*2+0] = axis[2]*axis[0]*icos - axis[1]*sine;
    rotmat[n*2+1] = axis[1]*axis[2]*icos + axis[0]*sine;
    rotmat[n*2+2] = axis[2]*axis[2]*icos + cosine;
}

/******************************************************************************/

void get_shiftmat_axis_angle(const double axis[3], const double angle,
        const bool forward, double* shiftmat){

    double axis_[3];
    double rotmat[9]; /* 3 x 3 matrix as 1-D array */

    for (int i=0; i<3; ++i){
        axis_[i] = -axis[i];
    }

    get_rotmat_axis_angle(axis_, angle, rotmat);

    if (!forward) {
        transpose(3, 3, rotmat, shiftmat);
    }
    else {
        memcpy(shiftmat, rotmat, 9*sizeof(double));
    }
}

/******************************************************************************/

/* Rotates vectors about axis by angle. v is a pointer to an array of
 * vectors, i.e. a n x 3 array. */
void rotate_vectors_axis_angle(const int n, const double* v, const double axis[3],
        const double angle, double* vrot){

    double rotmat[9]; /* 3 x 3 matrix as 1-D array */

    get_rotmat_axis_angle(axis, angle, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void shift_vectors_axis_angle (const int n, const double* v, 
        const double axis[3], const double angle, const bool forward,
        double* vshift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    get_shiftmat_axis_angle(axis, angle, forward, shiftmat);
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void shift_tensor2_axis_angle (const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    get_shiftmat_axis_angle(axis, angle, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void shift_tensor3_axis_angle (const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_axis_angle(axis, angle, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void axis_angle_to_dcm(const double axis[3], const double angle,
        double* dcm){

    get_shiftmat_axis_angle(axis, angle, true, dcm);
}

/******************************************************************************/

void axis_angle_to_euler(const double axis[3], const double angle, 
        enum EulangSeq seq, const bool world, double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_axis_angle(axis, angle, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
}

/******************************************************************************/

void axis_angle_to_quat(const double axis[3], const double angle, double q[4]){

    double half_angle = 0.5*angle;
    double sin_half_angle = sin(half_angle);
    double nrm;

    q[0] = cos(half_angle);
    q[1] = sin_half_angle*axis[0];
    q[2] = sin_half_angle*axis[1];
    q[3] = sin_half_angle*axis[2];

    /* Normalize the quaternion */
    nrm = norm(4, q);
    for (int i=0; i<4; ++i){
        q[i] /= nrm;
    }
}

/******************************************************************************/
                         /* DIRECTION COSINE MATRIX */
/******************************************************************************/


    /* Returns the direction cosine matrix of axes(i.e. frame) B w.r.t.      */
    /* axes(i.e. frame) A.                                                   */
    /*                                                                       */
    /* Parameters                                                            */
    /* ----------                                                            */
    /* A : 3 x 3 matrix                                                      */
    /*     The rows of A represent the orthonormal basis vectors of frame A. */
    /*                                                                       */
    /* B : 3 x 3 matrix                                                      */
    /*     The rows of B represent the orthonormal basis vectors of frame B. */
    /*                                                                       */
    /* Returns                                                               */
    /* -------                                                               */
    /* 3 x 3 matrix                                                          */
    /*     The dcm of frame B w.r.t. frame A.                                */


void dcm_from_axes(const double* A, const double* B, double* dcm) {

    double A_transpose[9];  /* 3 x 3 matrix as 1-D array */

    transpose(3, 3, A, A_transpose);
    mm_mult(3, 3, B, 3, A_transpose, dcm);
}

/******************************************************************************/

void get_rotmat_dcm(const double* dcm, double* rotmat){

    transpose(3, 3, dcm, rotmat);
}

/******************************************************************************/

void get_shiftmat_dcm(const double* dcm, const bool forward, double* shiftmat){

    if (!forward) {
        transpose(3, 3, dcm, shiftmat);
    }
    else {
        memcpy(shiftmat, dcm, 9*sizeof(double));
    }
}

/******************************************************************************/

void rotate_vectors_dcm(const int n, const double* v, const double* dcm,
        double* vrot){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_dcm(dcm, rotmat);
    rotate_vectors (n, v, rotmat, vrot);
}

/******************************************************************************/

void shift_vectors_dcm(const int n, const double* v, const double* dcm,
        const bool forward, double* vshift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_dcm(dcm, forward, shiftmat);
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void shift_tensor2_dcm(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_dcm(dcm, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void shift_tensor3_dcm(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_dcm(dcm, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void dcm_to_quat(const double* dcm, double q[4]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */
    double axis[3];
    double angle = 0.0;

    get_rotmat_dcm(dcm, rotmat);
    extract_axis_angle_from_rotmat(rotmat, axis, &angle);
    axis_angle_to_quat(axis, angle, q);
}

/******************************************************************************/

void dcm_to_euler(const double* dcm, enum EulangSeq seq, const bool world,
        double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_dcm(dcm, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
}

/******************************************************************************/

void dcm_to_axis_angle(const double* dcm, double axis[3], double* angle){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_dcm(dcm, rotmat);
    extract_axis_angle_from_rotmat(rotmat, axis, angle);
}

/******************************************************************************/
                               /* EULER ANGLES */
/******************************************************************************/

void get_rotmat_euler(const double euler[3], enum EulangSeq seq, const bool world,
        double* rotmat){

    rotmat_euler(euler, seq, world, rotmat);
}

/******************************************************************************/

void get_shiftmat_euler(const double euler[3], enum EulangSeq seq, const bool world,
        const bool forward, double* shiftmat){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    rotmat_euler(euler, seq, world, rotmat);

    if (forward){
        transpose(3, 3, rotmat, shiftmat);
    }
    else {
        memcpy(shiftmat, rotmat, 9*sizeof(double));
    }
}

/******************************************************************************/

void rotate_vectors_euler(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, double* vrot) {

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_euler(euler, seq, world, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void shift_vectors_euler(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* vshift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_euler(euler, seq, world, forward, shiftmat);
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void shift_tensor2_euler(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_euler(euler, seq, world, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void shift_tensor3_euler(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_euler(euler, seq, world, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void euler_to_axis_angle(const double euler[3], enum EulangSeq seq, const bool world,
        double axis[3], double* angle){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_euler(euler, seq, world, rotmat);
    extract_axis_angle_from_rotmat(rotmat, axis, angle);
}

/******************************************************************************/

void euler_to_dcm(const double euler[3], enum EulangSeq seq, const bool world,
        double* dcm){

    get_shiftmat_euler(euler, seq, world, true, dcm);
}

/******************************************************************************/

void euler_to_euler(const double euler[3], enum EulangSeq seq, const bool world,
        enum EulangSeq to_seq, const bool to_world, double to_euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_euler(euler, seq, world, rotmat);
    factor_rotmat(rotmat, to_seq, to_world, to_euler);
}

/******************************************************************************/

void euler_to_quat(const double euler[3], enum EulangSeq seq, const bool world,
        double q[4]){

    double axis[3];
    double angle;

    euler_to_axis_angle(euler, seq, world, axis, &angle);
    axis_angle_to_quat(axis, angle, q);
}

/******************************************************************************/
                               /* QUATERNIONS */
/******************************************************************************/

void get_rand_quat(double q[4]){

    double axis[3];
    double angle;

    get_rand_axis_angle(axis, &angle);
    axis_angle_to_quat(axis, angle, q);
}

/******************************************************************************/

void get_identity_quat(double q[4]){

    q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0;
}

/******************************************************************************/

/* Performs in-place conjugation of a quaternion */
void conjugate_quat(double q[4]){

    q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
}

/******************************************************************************/

/* Performs in-place inversion of a quaternion */
void invert_quat(double q[4]){

    conjugate_quat(q);
}

/******************************************************************************/

/* Performs in-place normalization of a quaternion */
void normalize_quat(double q[4]){

    double nrm = norm(4, q);

    for (int i=0; i<4; ++i){
        q[i] /= nrm;
    }
}

/******************************************************************************/

/* Returns true if a quaternion is normalized */
bool quat_is_normalized(const double q[4]){

    double nrm = norm(4, q);

    if (isclose(nrm, 1.0, 1e-14, 0.0)){
        return true;
    }
    else {
        return false;
    }
}

/******************************************************************************/

/* Calculates the product of two unit quaternions */
void get_quat_prod(const double p[4], const double q[4], double pq[4]){

    pq[0] = p[0]*q[0] - p[1]*q[1] -p [2]*q[2] -p [3]*q[3];
    pq[1] = p[1]*q[0] + p[0]*q[1] -p [3]*q[2] +p [2]*q[3];
    pq[2] = p[2]*q[0] + p[3]*q[1] +p [0]*q[2] -p [1]*q[3];
    pq[3] = p[3]*q[0] - p[2]*q[1] +p [1]*q[2] +p [0]*q[3];

    normalize_quat(pq);
}

/******************************************************************************/

double get_angle_between_quat(const double p[4], const double q[4]){

    return acos(dot(4, p, q));
}

/******************************************************************************/

void interpolate_quat(const double q1[4], const double q2[4], const double t,
        double q[4]){

    double theta = get_angle_between_quat(q1, q2);
    double sin_theta = sin(theta);
    double sin_t_theta = sin(t*theta);
    double sin_it_theta = sin((1.0-t)*theta);

    for (int i=0; i<4; ++i){
        q[i] = ( sin_it_theta*q1[i] + sin_t_theta*q2[i] )/sin_theta;
    }
}

/******************************************************************************/

void quat_deriv_to_ang_vel_mat(const double q[4], double* mat){

    /* mat is 3 x 4 */
    const int n = 4;

    mat[n*0+0] = -q[1]; mat[n*0+1] =  q[0]; mat[n*0+2] = -q[3]; mat[n*0+3] =  q[2];
    mat[n*1+0] = -q[2]; mat[n*1+1] =  q[3]; mat[n*1+2] =  q[0]; mat[n*1+3] = -q[1];
    mat[n*2+0] = -q[3]; mat[n*2+1] = -q[2]; mat[n*2+2] =  q[1]; mat[n*2+3] =  q[0];

    for (int i=0; i<3; ++i){
        for (int j=0; j<4; ++j){
            mat[4*i+j] *= 2.0;
        }
    }
}

/******************************************************************************/

void ang_vel_to_quat_deriv_mat(const double q[4], double* mat){

    /* mat is 4 x 3 */
    const int n = 3;

    mat[n*0+0] = -q[1]; mat[n*0+1] = -q[2]; mat[n*0+2] = -q[3];
    mat[n*1+0] =  q[0]; mat[n*1+1] =  q[3]; mat[n*1+2] = -q[2];
    mat[n*2+0] = -q[3]; mat[n*2+1] =  q[0]; mat[n*2+2] =  q[1];
    mat[n*3+0] =  q[2]; mat[n*3+1] = -q[1]; mat[n*3+2] =  q[0];

    for (int i=0; i<4; ++i){
        for (int j=0; j<3; ++j){
            mat[3*i+j] *= 0.5;
        }
    }
}

/******************************************************************************/

void quat_deriv_to_ang_vel(const double q[4], const double qdot[4], 
        double ang_vel[3]){

    double mat[12]; /* 3 x 4 matrix as 1-D array */

    quat_deriv_to_ang_vel_mat(q, mat);
    mv_mult(3, 4, mat, qdot, ang_vel);
}

/******************************************************************************/

void ang_vel_to_quat_deriv(const double q[4], const double ang_vel[3],
        double qdot[4]){

    double mat[12]; /* 4 x 3 matrix as 1-D array */

    ang_vel_to_quat_deriv_mat(q, mat);
    mv_mult(4, 3, mat, ang_vel, qdot);
}

/******************************************************************************/

void get_rotmat_quat(const double q[4], double* rotmat){

    const int n = 3;

    double q0sq = q[0]*q[0];
    double q1sq = q[1]*q[1];
    double q2sq = q[2]*q[2];
    double q3sq = q[3]*q[3];
    double q0q1 = q[0]*q[1];
    double q0q2 = q[0]*q[2];
    double q0q3 = q[0]*q[3];
    double q1q2 = q[1]*q[2];
    double q1q3 = q[1]*q[3];
    double q2q3 = q[2]*q[3];

    rotmat[n*0+0] = 2*(q0sq + q1sq) - 1.0;
    rotmat[n*0+1] = 2*(q1q2 - q0q3);
    rotmat[n*0+2] = 2*(q1q3 + q0q2);
    rotmat[n*1+0] = 2*(q1q2 + q0q3);
    rotmat[n*1+1] = 2*(q0sq + q2sq) - 1.0;
    rotmat[n*1+2] = 2*(q2q3 - q0q1);
    rotmat[n*2+0] = 2*(q1q3 - q0q2);
    rotmat[n*2+1] = 2*(q2q3 + q0q1);
    rotmat[n*2+2] = 2*(q0sq + q3sq) - 1.0;
}

/******************************************************************************/

void get_shiftmat_quat(const double q[4], const bool forward,
        double* shiftmat){

    double conj_q[4];

    if (forward){
        memcpy(conj_q, q, 4*sizeof(double));
        conjugate_quat(conj_q);
        get_rotmat_quat(conj_q, shiftmat);
    }
    else {
        get_rotmat_quat(q, shiftmat);
    }
}

/******************************************************************************/

void rotate_vectors_quat(const int n, const double* v, const double q[4],
        double* vrot){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_quat(q, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void shift_vectors_quat(const int n, const double* v, const double q[4],
        const bool forward, double* vshift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_quat(q, forward, shiftmat);
    shift_vectors(n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void shift_tensor2_quat(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_quat(q, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void shift_tensor3_quat(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    get_shiftmat_quat(q, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void quat_to_axis_angle(const double q[4], double axis[3], double* angle){

    double sine = sqrt(1.0-q[0]*q[0]);

    *angle = 2*acos(q[0]);

    if (*angle > 0.0) {
        if (*angle < M_PI){
            axis[0] = q[1]/sine;
            axis[1] = q[2]/sine;
            axis[2] = q[3]/sine;
        }
        else {
            double rotmat[9];  /* 3 x 3 matrix as 1-D array */
            get_rotmat_quat(q, rotmat);
            extract_axis_angle_from_rotmat(rotmat, axis, angle);
        }
    }
    else {
        axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
    }
    fix_axis_angle(axis, angle, true);
}

/******************************************************************************/

void quat_to_dcm(const double q[4], double* dcm){

    get_shiftmat_quat(q, true, dcm);
}

/******************************************************************************/

void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    get_rotmat_quat(q, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
}

/******************************************************************************/
                                  /* OTHERS */
/******************************************************************************/

/* In-place translations of a set of vectors */
void translate(const int n, double* v, const double delta[3]){

    for (int i=0; i<n; ++i){
        for (int j=0; j<3; ++j){
            v[3*i+j] += delta[j];
        }
    }
}

/******************************************************************************/

bool mat_is_rotmat(const double* mat){

    double identity[9];  /* 3 x 3 matrix as 1-D array */
    double mt[9];        /* 3 x 3 matrix as 1-D array */
    double mmt[9];       /* 3 x 3 matrix as 1-D array */

    double determinant = det(3, mat);
    bool det_is_one = isclose(determinant, 1.0, 1e-12, 1e-12);

    eye(3, identity);

    transpose(3, 3, mat, mt);
    mm_mult(3, 3, mat, 3, mt, mmt);

    bool is_orthogonal = allclose(9, mmt, identity, 1e-12, 1e-12);

    return (is_orthogonal && det_is_one);
}

/******************************************************************************/

bool mat_is_dcm(const double* mat){

    return mat_is_rotmat(mat);
}

/******************************************************************************/
