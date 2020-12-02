#ifndef  ROTATION_H
#define ROTATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "rot_shift_mat.h"
#include "rot_eulang.h"
#include <stdbool.h>


/******************************************************************************/
                                /* AXIS ANGLE */
/******************************************************************************/

void fix_axis_angle(double axis[3], double* angle, const bool normalize);

void get_rand_axis_angle(double axis[3], double* angle);

void get_rotmat_axis_angle(const double axis[3], const double angle,
        double* rotmat);

void get_shiftmat_axis_angle(const double axis[3], const double angle,
        const bool forward, double* shiftmat);

void rotate_vectors_axis_angle(const int n, const double* v, const double axis[3],
        const double angle, double* vrot);

void shift_vectors_axis_angle (const int n, const double* v, 
        const double axis[3], const double angle, const bool forward,
        double* vshift);

void shift_tensor2_axis_angle (const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

void shift_tensor3_axis_angle (const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift);

void axis_angle_to_dcm(const double axis[3], const double angle,
        double* dcm);

void axis_angle_to_euler(const double axis[3], const double angle, 
        enum EulangSeq seq, const bool world, double euler[3]);

void axis_angle_to_quat(const double axis[3], const double angle, double q[4]);


/******************************************************************************/
                         /* DIRECTION COSINE MATRIX */
/******************************************************************************/

void dcm_from_axes(const double* A, const double* B, double* dcm);

void get_rotmat_dcm(const double* dcm, double* rotmat);

void get_shiftmat_dcm(const double* dcm, const bool forward, double* shiftmat);

void rotate_vectors_dcm(const int n, const double* v, const double* dcm,
        double* vrot);

void shift_vectors_dcm(const int n, const double* v, const double* dcm,
        const bool forward, double* vshift);

void shift_tensor2_dcm(const double* A, const double* dcm, const bool forward,
        double* Ashift);

void shift_tensor3_dcm(const double* A, const double* dcm, const bool forward,
        double* Ashift);

void dcm_to_quat(const double* dcm, double q[4]);

void dcm_to_euler(const double* dcm, enum EulangSeq seq, const bool world,
        double euler[3]);

void dcm_to_axis_angle(const double* dcm, double axis[3], double* angle);


/******************************************************************************/
                               /* EULER ANGLES */
/******************************************************************************/

void get_rotmat_euler(const double euler[3], enum EulangSeq seq, const bool world,
        double* rotmat);

void get_shiftmat_euler(const double euler[3], enum EulangSeq seq, const bool world,
        const bool forward, double* shiftmat);

void rotate_vectors_euler(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, double* vrot);

void shift_vectors_euler(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* vshift);

void shift_tensor2_euler(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

void shift_tensor3_euler(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward, double* Ashift);

void euler_to_axis_angle(const double euler[3], enum EulangSeq seq, const bool world,
        double axis[3], double* angle);

void euler_to_dcm(const double euler[3], enum EulangSeq seq, const bool world,
        double* dcm);

void euler_to_euler(const double euler[3], enum EulangSeq seq, const bool world,
        enum EulangSeq to_seq, const bool to_world, double to_euler[3]);

void euler_to_quat(const double euler[3], enum EulangSeq seq, const bool world,
        double q[4]);

/******************************************************************************/
                               /* QUATERNIONS */
/******************************************************************************/

void get_rand_quat(double q[4]);

void get_identity_quat(double q[4]);

/* Performs in-place conjugation of a quaternion */
void conjugate_quat(double q[4]);

/* Performs in-place inversion of a quaternion */
void invert_quat(double q[4]);

/* Performs in-place normalization of a quaternion */
void normalize_quat(double q[4]);

/* Returns true if a quaternion is normalized */
bool quat_is_normalized(const double q[4]);

/* Calculates the product of two unit quaternions */
void get_quat_prod(const double p[4], const double q[4], double pq[4]);

double get_angle_between_quat(const double p[4], const double q[4]);

void interpolate_quat(const double q1[4], const double q2[4], const double t,
        double q[4]);

/* mat is 3 x 4 */
void quat_deriv_to_ang_vel_mat(const double q[4], double* mat);

/* mat is 4 x 3 */
void ang_vel_to_quat_deriv_mat(const double q[4], double* mat);

void quat_deriv_to_ang_vel(const double q[4], const double qdot[4], 
        double ang_vel[3]);

void ang_vel_to_quat_deriv(const double q[4], const double ang_vel[3],
        double qdot[4]);

void get_rotmat_quat(const double q[4], double* rotmat);

void get_shiftmat_quat(const double q[4], const bool forward,
        double* shiftmat);

void rotate_vectors_quat(const int n, const double* v, const double q[4],
        double* vrot);

void shift_vectors_quat(const int n, const double* v, const double q[4],
        const bool forward, double* vshift);

void shift_tensor2_quat(const double* A, const double q[4],
        const bool forward, double* Ashift);

void shift_tensor3_quat(const double* A, const double q[4],
        const bool forward, double* Ashift);

void quat_to_axis_angle(const double q[4], double axis[3], double* angle);

void quat_to_dcm(const double q[4], double* dcm);

void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        double euler[3]);


/******************************************************************************/
                                  /* OTHERS */
/******************************************************************************/

/* In-place translations of a set of vectors */
void translate(const int n, double* v, const double delta[3]);

bool mat_is_rotmat(const double* mat);

bool mat_is_dcm(const double* mat);

#ifdef __cplusplus
}
#endif

#endif /* ROTATION_H */
