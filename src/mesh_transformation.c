#include "mesh_transformation.h"
#include "ran_num.h"
#include "utils_math.h"
#include <stdlib.h>

/******************************************************************************/

/* Translate a mesh by a vector dx */
void mesh_translate(Mesh* msh, const double* dx){

    int nc;
    int num_vertices = mesh_get_num_vertices(msh);

    for (int ivert=0; ivert < num_vertices; ++ivert){
        double* coords = mesh_get_coords(msh, ivert, &nc);
        for (int j=0; j < 3; ++j){
            coords[j] += dx[j];
        }
    }
}

/******************************************************************************/

/* Scale mesh by a given scaling factors */
void mesh_scale(Mesh* msh, const double sx, const double sy, const double sz)
{
    int num_vertices = mesh_get_num_vertices(msh);
    double* coordinates = mesh_get_coordinates(msh, NULL);

    for (int ivert=0; ivert < num_vertices; ++ivert){
        coordinates[3*ivert+0] *= sx;
        coordinates[3*ivert+1] *= sy;
        coordinates[3*ivert+2] *= sz;
    }
}

/******************************************************************************/

/* Rotate a mesh by an angle about a given axis*/
void mesh_rotate_axis_angle(Mesh* msh, const double* axis, const double angle){

    int nc;

    int num_vertices = mesh_get_num_vertices(msh);

    double* coordinates = mesh_get_coordinates(msh, &nc);

    double* crot = malloc(3*num_vertices*sizeof(double));

    rotate_vectors_axis_angle(num_vertices, coordinates, axis, angle, crot);

    mesh_set_coordinates(msh, crot);

    free(crot);

}

/******************************************************************************/

/* Rotate a mesh by euler angle */
void mesh_rotate_eulang(Mesh* msh, const double* euler, enum EulangSeq seq,
        const bool world){

    int nc;

    int num_vertices = mesh_get_num_vertices(msh);

    double* coordinates = mesh_get_coordinates(msh, &nc);

    double* crot = malloc(3*num_vertices*sizeof(double));

    rotate_vectors_euler(num_vertices, coordinates, euler, seq, world, crot);

    mesh_set_coordinates(msh, crot);

    free(crot);

}

/******************************************************************************/

void mesh_perturb(Mesh* msh, const double dir[3], const double amplitude)
{
    double lb = -0.5*amplitude;
    double ub =  0.5*amplitude;
    int num_coordinates;
    double* coordinates = mesh_get_coordinates(msh, &num_coordinates);
    int num_vertices = num_coordinates/3;

    double* prtrb = malloc(num_coordinates*sizeof(double));

    RandomStream stream = init_stream(0);

    duni(stream, lb, ub, num_coordinates, prtrb);

    for (int ivert=0; ivert < num_vertices; ++ivert){
        for (int j=0; j < 3; ++j){
            coordinates[3*ivert+j] += dir[j]*prtrb[3*ivert+j];
        }
    }

    delete_stream(&stream);

    free(prtrb);
}

/******************************************************************************/

void mesh_get_com(const Mesh* msh, double* com)
{
    int num_coordinates;
    int num_verts = mesh_get_num_vertices(msh);
    double* coordinates = mesh_get_coordinates(msh, &num_coordinates);

    zero_out(3, com);

    for (int ivert=0; ivert < num_verts; ++ivert){
        for (int j=0; j < 3; ++j){
            com[j] += coordinates[3*ivert+j];
        }
    }
    for (int j=0; j < 3; ++j){
        com[j] /= num_verts;
    }
}

/******************************************************************************/

void mesh_pull_back(Mesh* msh)
{
    double com[3];

    mesh_get_com(msh, com);

    for (int j=0; j < 3; ++j){
        com[j] = -com[j];
    }

    mesh_translate(msh, com);
}

/******************************************************************************/
