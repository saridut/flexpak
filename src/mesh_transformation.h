#ifndef MESH_TRANSFORMATION_H
#define MESH_TRANSFORMATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "mesh.h"
#include "rotation.h"

/* Translate a mesh by a vector dx */
void mesh_translate(Mesh* msh, const double* dx);

/* Rescale mesh by given scaling factors */
void mesh_scale(Mesh* msh, const double sx, const double sy, const double sz);

/* Rotate a mesh by an angle about a given axis*/
void mesh_rotate_axis_angle(Mesh* msh, const double* axis, const double angle);

/* Rotate a mesh by euler angle */
void mesh_rotate_eulang(Mesh* msh, const double* euler, enum EulangSeq seq,
        const bool world);

/*
 *  Add random perturbations to the vertex coordinates
 */
void mesh_perturb(Mesh* msh, const double dir[3], const double amplitude);


void mesh_get_com(const Mesh* msh, double* com);


void mesh_pull_back(Mesh* msh);


#ifdef __cplusplus
}
#endif

#endif //MESH_TRANSFORMATION_H
