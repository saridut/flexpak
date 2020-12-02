#ifndef MESH_H
#define MESH_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mesh_connectivity.h"
#include "mesh_geometry.h"
#include <stdbool.h>

struct Mesh{
    int tdim;
    char cell_type[8];
    MeshConnectivity** connectivities;
    MeshGeometry geometry;
    int* num_entities;
};

typedef struct Mesh Mesh;

/* Use later if needed:                                   */
/* enum entity {VERTEX=0, EDGE=1, FACE=2, CELL=3, ALL=4}; */

/* Creates an empty mesh of topological dimension tdim*/
Mesh mesh_create(const int tdim, const char* cell_type);


/* Copies a mesh */
void mesh_copy(Mesh* this, const Mesh* other);


/* Deletes a mesh */
void mesh_delete(Mesh* this);


/* Clears a mesh */
void mesh_clear(Mesh* this);


/* Adds geometry from a file */
void mesh_add_geometry_from_file(Mesh* this, const char* fn);


/* Adds geometry from an array */
void mesh_add_geometry_from_array(Mesh* this, const int num_vertices,
        const double* coordinates);


/* Adds connectivity from a file */
void mesh_add_connectivity_from_file(Mesh* this, const int d0, const int d1,
        const char* fn);


/* Gets the number of vertices */
int mesh_get_num_vertices(const Mesh* this);


/* Gets the total number of coordinates */
int mesh_get_num_coordinates(const Mesh* this);


/* Gets a pointer to all coordinates */
double* mesh_get_coordinates(const Mesh* this, int* num_coordinates);


/* Gets a pointer to the coordinates of vertex ivert */
double* mesh_get_coords(const Mesh* this, const int ivert, int* nc);


/* Gets a pointer to the coordinate i of vertex ivert; 0<=i<=2 */
double* mesh_get_coord(const Mesh* this, const int ivert, const int i);


/* Sets all coordinates */
void mesh_set_coordinates(Mesh* this, const double* coordinates);


/* Sets the coordinates of vertex ivert */
void mesh_set_coords(Mesh* this, const int ivert, const double* coords);


/* Sets coordinate i of vertex ivert; 0<=i<=2 */
void mesh_set_coord(Mesh* this, const int ivert, const int i, const double* coord);


/* Gets the number of mesh entities of topological dimension d*/
int mesh_get_num_entities(const Mesh* this, const int d);


/* Gets the number of connections in connectivity (d0->d1) */
int mesh_get_num_connections(const Mesh* this, const int d0, const int d1);


/* Gets the number of connections for entity ient in connectivity (d0->d1) */
int mesh_get_num_conns(const Mesh* this, const int d0, const int d1,
        const int ient);


/* Gets a pointer to mesh connectivity for (d0->d1)  */
MeshConnectivity* mesh_get_connectivity(const Mesh* this, const int d0,
        const int d1);


/* Gets a pointer to connections for entity ient in (d0->d1) */
int* mesh_get_conns(const Mesh* this, const int d0,
        const int d1, const int ient, int* num_conns);


/* Gets a pointer to connection i for entity ient in (d0->d1) */
int* mesh_get_conn(const Mesh* this, const int d0,
        const int d1, const int ient, const int i);

/* Get coordinates for an entity ient of topological dimension d.    */
/* The coordinates are copied to a buffer with pointer coords.       */
/* If num_coords is not NULL, the num_coords contains the number     */
/* of coordinates. Note that num_coords is the total number of x, y, */
/* and z-coordinates -- it can by divided by three to loop over      */
/* incident vertices.                                                */
void mesh_get_ent_coords(const Mesh* this, const int d, const int ient,
        double* coords, int* num_coords);


/* Returns true if entity jent of topological dimension d1 is incident */
/* to entity ient of topological dimension d0.                         */
bool mesh_is_incident(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent);

/* Returns the local index of entity (d1, jent) within entity */
/* (d0, ient), where d0 > d1                                  */
int mesh_get_local_index(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent);


/* Returns the global index of entity (d1, ient:jent) within */
/* entity (d0, ient), where d0 > d1.                         */
int mesh_get_global_index(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent);


/* Writes mesh to an obj file */
void mesh_to_obj(const Mesh* this, const char* fn);

/* Writes mesh to a vtk file */
void mesh_to_vtk(const Mesh* this, const char* fn);

/* Write to hdf5 buffer? */

/* void mesh_to_hdf5(Mesh* this, ...); */

#ifdef __cplusplus
}
#endif

#endif //MESH_H
