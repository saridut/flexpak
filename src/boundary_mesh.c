#include "boundary_mesh.h"
#include "dynamic_array.h"
#include <string.h>

/*****************************************************************************/

BoundaryMesh bmesh_create(){

    BoundaryMesh bmsh;
    
    bmsh.tdim = -1;
    bmsh.msh = NULL;
    bmsh.num_cells = 0;
    bmsh.num_vertices = 0;
    bmsh.cell_map = NULL;
    bmsh.vertex_map = NULL;

    return bmsh;
}

/*****************************************************************************/

void bmesh_init(BoundaryMesh* this, Mesh* msh){

    this->tdim = msh->tdim - 1;
    this->msh = msh;

    /* Total number of facets in the full mesh */
    int num_facets_msh = mesh_get_num_entities(msh, this->tdim);

    /* Loop over all facets of the full mesh to find boundary cells */
    iDynamicArray cm = idyar_create(0);
    for (int ifacet=0; ifacet < num_facets_msh; ++ifacet){
        int nc = mesh_get_num_conns(msh, this->tdim, msh->tdim, ifacet);
        if (nc == 1){
            idyar_append(&cm, ifacet);
        }
    }
    this->num_cells = idyar_size(&cm);

    this->cell_map = malloc(this->num_cells*sizeof(int));
    memcpy(this->cell_map, idyar_data(&cm), this->num_cells*sizeof(int));
    idyar_delete(&cm);

    /* Loop over boundary facets to collect boundary vertices */
    iDynamicArray vm = idyar_create(0);
    int nc;
    for (int i=0; i < this->num_cells; ++i){
        int i_msh = this->cell_map[i];
        int* conns = mesh_get_conns(msh, this->tdim, 0, i_msh, &nc);
        for (int j=0; j < nc; ++j){
            idyar_append(&vm, conns[j]);
        }
    }
    idyar_unique(&vm);

    this->num_vertices = idyar_size(&vm);
    this->vertex_map = malloc(this->num_vertices*sizeof(int));
    memcpy(this->vertex_map, idyar_data(&vm), this->num_vertices*sizeof(int));
    idyar_delete(&vm);
}

/*****************************************************************************/

/* Deletes a boundary mesh */
void bmesh_delete(BoundaryMesh* this){

    this->msh = NULL;
    free(this->cell_map);
    free(this->vertex_map);
}

/*****************************************************************************/

/* Returns the number of boundary cells */
int bmesh_get_num_cells(const BoundaryMesh* this){

    return this->num_cells;
}

/*****************************************************************************/

/* Returns the number of boundary vertices */
int bmesh_get_num_vertices(const BoundaryMesh* this){

    return this->num_vertices;
}

/*****************************************************************************/

/* Gets a pointer to the coordinates of boundary vertex ivert */
double* bmesh_get_coords(const BoundaryMesh* this, const int ivert, int* nc){

    return mesh_get_coords(this->msh, this->vertex_map[ivert], nc);
}

/*****************************************************************************/

/* Gets a pointer to the coordinate i of boundary vertex ivert; 0<=i<=2 */
double* bmesh_get_coord(const BoundaryMesh* this, const int ivert, 
        const int i){

    return mesh_get_coord(this->msh, this->vertex_map[ivert], i);
}

/*****************************************************************************/

/* Sets the coordinates of boundary vertex ivert */
void bmesh_set_coords(BoundaryMesh* this, const int ivert,
        const double* coords) {

    mesh_set_coords(this->msh, this->vertex_map[ivert], coords);

}

/*****************************************************************************/

/* Sets coordinate i of boundary vertex ivert; 0<=i<=2 */
void bmesh_set_coord(BoundaryMesh* this, const int ivert, const int i,
        const double* coord){

    mesh_set_coord(this->msh, this->vertex_map[ivert], i, coord);
}

/*****************************************************************************/

/* Gets the number of connections for entity ient of boundary mesh in */
/* connectivity (d0->d1)                                              */
int bmesh_get_num_conns(const BoundaryMesh* this, const int d0, const int d1,
        const int ient){

    return mesh_get_num_conns(this->msh, d0, d1, this->cell_map[ient]);
}

/*****************************************************************************/

/* Gets a pointer to connections for entity ient of boundary mesh in */
/* connectivity (d0->d1)                                             */
int* bmesh_get_conns(const BoundaryMesh* this, const int d0,
        const int d1, const int ient, int* num_conns){

    return mesh_get_conns(this->msh, d0, d1, this->cell_map[ient], num_conns);
}

/*****************************************************************************/

/* Gets a pointer to connection i for entity ient of boundary mesh in (d0->d1) */
int* bmesh_get_conn(const BoundaryMesh* this, const int d0,
        const int d1, const int ient, const int i){

    return mesh_get_conn(this->msh, d0, d1, this->cell_map[ient], i);
}

/*****************************************************************************/

/* Returns the corresponding entity of topological dimension d in the full mesh */
int bmesh_get_mapped_entity(const BoundaryMesh* this, const int d,
        const int ient){

    if (d == 0){
        return this->vertex_map[ient];
    }
    else if (d == 1){
        return this->cell_map[ient];
    }
    else {
        return -1;
    }
}

/*****************************************************************************/
