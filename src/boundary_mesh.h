#ifndef BOUNDARY_MESH_H
#define BOUNDARY_MESH_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mesh.h"

struct BoundaryMesh{
    int tdim;
    Mesh* msh;
    int num_cells;
    int num_vertices;
    int* cell_map;
    int* vertex_map;
};

typedef struct BoundaryMesh BoundaryMesh;


/* Creates an empty boundary mesh */
BoundaryMesh bmesh_create();


/* Initializes from a full mesh */
void bmesh_init(BoundaryMesh* this, Mesh* msh);


/* Deletes a boundary mesh */
void bmesh_delete(BoundaryMesh* this);


/* Returns the number of boundary cells */
int bmesh_get_num_cells(const BoundaryMesh* this);


/* Returns the number of boundary vertices */
int bmesh_get_num_vertices(const BoundaryMesh* this);


/* Gets a pointer to the coordinates of boundary vertex ivert */
double* bmesh_get_coords(const BoundaryMesh* this, const int ivert, int* nc);


/* Gets a pointer to the coordinate i of boundary vertex ivert; 0<=i<=2 */
double* bmesh_get_coord(const BoundaryMesh* this, const int ivert, const int i);


/* Sets the coordinates of boundary vertex ivert */
void bmesh_set_coords(BoundaryMesh* this, const int ivert, const double* coords);


/* Sets coordinate i of boundary vertex ivert; 0<=i<=2 */
void bmesh_set_coord(BoundaryMesh* this, const int ivert, const int i,
        const double* coord);


/* Gets the number of connections for entity ient of boundary mesh in */
/* connectivity (d0->d1)                                              */
int bmesh_get_num_conns(const BoundaryMesh* this, const int d0, const int d1,
        const int ient);


/* Gets a pointer to connections for entity ient of boundary mesh in */
/* connectivity (d0->d1)                                             */
int* bmesh_get_conns(const BoundaryMesh* this, const int d0,
        const int d1, const int ient, int* num_conns);


/* Gets a pointer to connection i for entity ient of boundary mesh in (d0->d1) */
int* bmesh_get_conn(const BoundaryMesh* this, const int d0,
        const int d1, const int ient, const int i);


/* Returns the corresponding entity of topological dimension d in the full mesh */
int bmesh_get_mapped_entity(const BoundaryMesh* this, const int d, const int ient);


/* Write cell_map to file (Not implemented) */
//void bmesh_cmap_to_file(const BoundaryMesh* this);

#ifdef __cplusplus
}
#endif

#endif //BOUNDARY_MESH_H
