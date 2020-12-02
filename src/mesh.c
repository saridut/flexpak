#include "mesh.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/******************************************************************************/

/* Creates an empty mesh of topological dimension tdim*/
Mesh mesh_create(const int tdim, const char* cell_type){
    Mesh msh;

    msh.tdim = tdim;
    assert(strlen(cell_type) < 8); /* See struct definition */
    strcpy(msh.cell_type, cell_type);

    msh.connectivities = malloc((tdim+1)*sizeof(MeshConnectivity*));
    for (int i=0; i < (tdim+1); ++i){
        msh.connectivities[i] = malloc((tdim+1)*sizeof(MeshConnectivity));
        for (int j=0; j < (tdim+1); ++j){
            msh.connectivities[i][j] = meshcon_create(i, j);
        }
    }

    msh.geometry = meshgeom_create();
    msh.num_entities = calloc(tdim+1, sizeof(int));

    return msh;
}

/******************************************************************************/

/* Copies a mesh of same topological dimension and cell_type */
void mesh_copy(Mesh* this, const Mesh* other){

    assert(this->tdim == other->tdim);
    assert(strcmp(this->cell_type, other->cell_type)==0);

    for (int i=0; i < (this->tdim+1); ++i){
        for (int j=0; j < (this->tdim+1); ++j){
            meshcon_copy(&this->connectivities[i][j],
                    &other->connectivities[i][j]);
        }
    }

    meshgeom_copy(&this->geometry, &other->geometry);
    memcpy(this->num_entities, other->num_entities, (this->tdim+1)*sizeof(int));
}

/******************************************************************************/

/* Deletes a mesh */
void mesh_delete(Mesh* this){

    /* Delete all mesh connectivities */
    for (int i=0; i < (this->tdim+1); ++i){
        for (int j=0; j < (this->tdim+1); ++j){
            meshcon_delete(&this->connectivities[i][j]);
        }
    }

    /* Free the connectivities */
    for (int i=0; i < (this->tdim+1); ++i){
        free(this->connectivities[i]);
    }
    free(this->connectivities);

    meshgeom_delete(&this->geometry);
    free(this->num_entities);
}

/******************************************************************************/

/* Clears a mesh (deletes mesh data for connectivities and geometry) */
void mesh_clear(Mesh* this){

    for (int i=0; i < (this->tdim+1); ++i){
        for (int j=0; j < (this->tdim+1); ++j){
            meshcon_clear(&this->connectivities[i][j]);
        }
    }

    meshgeom_clear(&this->geometry);

    for (int i=0; i < (this->tdim+1); ++i){
        this->num_entities[i] = 0;
    }
}

/******************************************************************************/

/* Adds geometry from a file */
void mesh_add_geometry_from_file(Mesh* this, const char* fn){

    meshgeom_from_file(&this->geometry, fn);
}

/******************************************************************************/

/* Adds geometry from an array */
void mesh_add_geometry_from_array(Mesh* this, const int num_vertices,
        const double* coordinates){

    meshgeom_from_array(&this->geometry, num_vertices, coordinates);
}

/******************************************************************************/

/* Adds connectivity from a file */
void mesh_add_connectivity_from_file(Mesh* this, const int d0, const int d1,
        const char* fn){

    meshcon_from_file(&this->connectivities[d0][d1], fn);

    int ne = meshcon_get_num_entities(&this->connectivities[d0][d1]);

    if (this->num_entities[d0] == 0){
        this->num_entities[d0] = ne ;
    }
    else {
        assert (this->num_entities[d0] == ne);
    }
}

/******************************************************************************/

/* Gets the number of vertices */
int mesh_get_num_vertices(const Mesh* this){

    return meshgeom_get_num_vertices(&this->geometry);
}

/******************************************************************************/

/* Gets the total number of coordinates */
int mesh_get_num_coordinates(const Mesh* this){

    return meshgeom_get_num_coordinates(&this->geometry);
}

/******************************************************************************/

/* Gets a pointer to all coordinates */
double* mesh_get_coordinates(const Mesh* this, int* num_coordinates){

    return meshgeom_get_coordinates(&this->geometry, num_coordinates);
}

/******************************************************************************/

/* Gets a pointer to the coordinates of vertex ivert */
double* mesh_get_coords(const Mesh* this, const int ivert, int* nc){

     return meshgeom_get_coords(&this->geometry, ivert, nc);
}

/******************************************************************************/

/* Gets a pointer to the coordinate i of vertex ivert; 0<=i<=2 */
double* mesh_get_coord(const Mesh* this, const int ivert, const int i){

    return meshgeom_get_coord(&this->geometry, ivert, i);
}

/******************************************************************************/

/* Sets all coordinates */
void mesh_set_coordinates(Mesh* this, const double* coordinates){

    meshgeom_set_coordinates(&this->geometry, coordinates);
}

/******************************************************************************/

/* Sets the coordinates of vertex ivert */
void mesh_set_coords(Mesh* this, const int ivert, const double* coords){

    meshgeom_set_coords(&this->geometry, ivert, coords);
}

/******************************************************************************/

/* Sets coordinate i of vertex ivert; 0<=i<=2 */
void mesh_set_coord(Mesh* this, const int ivert, const int i,
        const double* coord){

    meshgeom_set_coord(&this->geometry, ivert, i, coord);
}

/******************************************************************************/

/* Get number of mesh entities of topological dimension d*/
int mesh_get_num_entities(const Mesh* this, const int d){

    return this->num_entities[d];
}

/******************************************************************************/

/* Gets the number of connections in connectivity (d0->d1) */
int mesh_get_num_connections(const Mesh* this, const int d0, const int d1){

    return meshcon_get_num_connections(&this->connectivities[d0][d1]);
}

/******************************************************************************/

/* Gets the number of connections for entity ient in connectivity (d0->d1) */
int mesh_get_num_conns(const Mesh* this, const int d0, const int d1,
        const int ient){

    return meshcon_get_num_conns(&this->connectivities[d0][d1], ient);
}

/******************************************************************************/

/* Gets a pointer to mesh connectivity (d0->d1)  */
MeshConnectivity* mesh_get_connectivity(const Mesh* this, const int d0,
        const int d1){

    return &this->connectivities[d0][d1];
}

/******************************************************************************/

/* Gets a pointer to connections for entity ient in (d0->d1) */
int* mesh_get_conns(const Mesh* this, const int d0,
        const int d1, const int ient, int* num_conns){

    return meshcon_get_conns(&this->connectivities[d0][d1], ient, num_conns);
}

/******************************************************************************/

/* Gets a pointer to connection i for entity ient in (d0->d1) */
int* mesh_get_conn(const Mesh* this, const int d0,
        const int d1, const int ient, const int i){

    return meshcon_get_conn(&this->connectivities[d0][d1], ient, i);
}

/******************************************************************************/

/* Get coordinates for an entity ient of topological dimension d.    */
/* The coordinates are copied to a buffer with pointer coords.       */
/* If num_coords is not NULL, the num_coords contains the number     */
/* of coordinates. Note that num_coords is the total number of x, y, */
/* and z-coordinates -- it can by divided by three to loop over      */
/* incident vertices.                                                */
void mesh_get_ent_coords(const Mesh* this, const int d, const int ient,
        double* coords, int* num_coords){
    
    int nc;
    int* conns = mesh_get_conns(this, d, 0, ient, &nc);
    double* vert_ptr;

    for (int i=0; i < nc; ++i){
        vert_ptr = mesh_get_coords(this, conns[i], NULL);
        memcpy(coords+3*i, vert_ptr, 3*sizeof(double));
    }
    if (num_coords != NULL){
        *num_coords = 3*nc;
    }
}

/******************************************************************************/

/* Returns true if entity jent of topological dimension d1 is incident */
/* to entity ient of topological dimension d0.                         */
bool mesh_is_incident(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent){

    bool is_incident = false;
    int nc_ient;
    int* conns_ient;

    if (d0 < d1){
        return mesh_is_incident(this, d1, jent, d0, ient);
    }
    else if (d0 > d1){
        conns_ient = mesh_get_conns(this, d0, d1, ient, &nc_ient);
        for (int i=0; i < nc_ient; ++i){
            if (conns_ient[i] == jent){
                is_incident = true;
                break;
            }
        }
        return is_incident;
    }
    else { /* d0 == d1 */
        int nc_jent;
        int* conns_jent;
        if ( (d0==0) && (d1==0) ){
            /* Do the vertices belong to a common cell? */
            conns_ient = mesh_get_conns(this, 0, this->tdim, ient, &nc_ient);
            conns_jent = mesh_get_conns(this, 0, this->tdim, jent, &nc_jent);
        }
        else if ( (d0 > 0) && (d1 > 0) ){
            /* Are the entities incident to a common vertex? */
            conns_ient = mesh_get_conns(this, d0, 0, ient, &nc_ient);
            conns_jent = mesh_get_conns(this, d1, 0, jent, &nc_jent);
        }
        for (int i=0; i < nc_ient; ++i){
            for (int j=0; j < nc_jent; ++j){
                if (conns_ient[i] == conns_jent[j]){
                    is_incident = true;
                    goto jump;
                }
            }
        }
        jump: return is_incident;
    }
}

/******************************************************************************/

/* Returns the local index of entity (d1, jent) within entity */
/* (d0, ient), where d0 > d1                                  */
int mesh_get_local_index(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent)
{
    int nc;
    int* conns = mesh_get_conns(this, d0, d1, ient, &nc);
    int indx;
    for (int i=0; i < nc; ++i){
        if (conns[i] == jent){
            indx = i;
            break;
        }
    }
    return indx;
}

/******************************************************************************/

/* Returns the global index of entity (d1, ient:jent) within */
/* entity (d0, ient), where d0 > d1.                         */
int mesh_get_global_index(const Mesh* this, const int d0, const int ient,
        const int d1, const int jent){

    int nc;
    int* conns = mesh_get_conns(this, d0, d1, ient, &nc);

    return conns[jent];
}

/******************************************************************************/

/* Writes mesh to an obj file */
void mesh_to_obj(const Mesh* this, const char* fn){

    FILE* fp = fopen(fn, "w");

    int nverts = mesh_get_num_vertices(this);
    int nfacets = mesh_get_num_entities(this, 2); /* only triangles */
    int nc;

    /* Write vertices */
    for (int ivert=0; ivert < nverts; ++ivert){
        double* v = mesh_get_coords(this, ivert, &nc);
        fprintf (fp, "v");
        for (int j=0; j < nc; ++j){
            fprintf (fp, "  %f", v[j]);
        }
        fprintf (fp, "\n");
    }

    /* Write faces */
    for (int ifacet=0; ifacet < nfacets; ++ifacet){
        int* f = mesh_get_conns(this, this->tdim, 0, ifacet, &nc);
        fprintf (fp, "f");
        for (int j=0; j < nc; ++j){
            fprintf (fp, "  %d", f[j]+1); /* obj files begin counting from 1 */
        }
        fprintf (fp, "\n");
    }

    fclose(fp);
}

/******************************************************************************/

/* Writes mesh to a vtk file */
void mesh_to_vtk(const Mesh* this, const char* fn){

    FILE* fp = fopen(fn, "w");
    
    //Header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Sheet Mesh\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    int nverts = mesh_get_num_vertices(this);
    int nfacets = mesh_get_num_entities(this, 2); /* only triangles */
    int num_connections = mesh_get_num_connections(this, 2, 0);
    int nc;

    /* Write vertices */
    fprintf(fp, "POINTS %d double\n", nverts);
    for (int ivert=0; ivert < nverts; ++ivert){
        double* v = mesh_get_coords(this, ivert, &nc);
        for (int j=0; j < nc; ++j){
            fprintf (fp, "  %f", v[j]);
        }
        fprintf (fp, "\n");
    }

    /* Write faces */
    fprintf(fp, "CELLS %d  %d\n", nfacets, nfacets+num_connections);
    for (int ifacet=0; ifacet < nfacets; ++ifacet){
        int* f = mesh_get_conns(this, this->tdim, 0, ifacet, &nc);
        fprintf (fp, "%d", nc);
        for (int j=0; j < nc; ++j){
            fprintf (fp, "  %d", f[j]);
        }
        fprintf (fp, "\n");
    }

    fprintf(fp, "CELL_TYPES %d\n", nfacets);
    for (int ifacet=0; ifacet < nfacets; ++ifacet){
        fprintf(fp, "5\n"); //triangle is VTK type 5
    }

    fclose(fp);
}

/******************************************************************************/
