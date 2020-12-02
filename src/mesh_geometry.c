#include "mesh_geometry.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************/

/*!
 *  \brief Creates an empty MeshGeometry object
 */
MeshGeometry meshgeom_create(){

    MeshGeometry    mg;

    mg.num_vertices = 0;
    mg.num_coordinates = 0;
    mg.coordinates = NULL;
    mg.ind_pos = NULL;

    return mg;
}

/******************************************************************************/

/*!
 *  \brief Initializes from an array of coordinates.
 *
 * @param[in, out] this pointer to object
 * @param[in] num_vertices number of vertices
 * @param[in] coordinates  pointer to coordinates array
 */
void meshgeom_from_array(MeshGeometry* this, const int num_vertices,
        const double* coordinates){

    int num_coordinates = 3*num_vertices;

    this->num_vertices = num_vertices;
    this->num_coordinates = num_coordinates;

    this->ind_pos = malloc((num_vertices+1)*sizeof(int));
    this->coordinates = malloc(num_coordinates*sizeof(double));

    this->ind_pos[0] = 0;
    for (int i=0; i < num_vertices; ++i){
        this->ind_pos[i+1] = this->ind_pos[i] + 3;
    }

    memcpy(this->coordinates, coordinates, num_coordinates*sizeof(double));
}

/******************************************************************************/

/* Initializes from file */
void meshgeom_from_file(MeshGeometry* this, const char* fn){

    FILE* fp = fopen(fn, "r");
    int ivert;

    /* Read the header line */
    fscanf(fp, "%d %d", &this->num_vertices, &this->num_coordinates);
    
    /* Allocate memory for arrays ind_pos and coordinates */
    this->ind_pos = malloc((this->num_vertices+1)*sizeof(int));
    this->coordinates = malloc(this->num_coordinates*sizeof(double));

    /* Loop over all lines. For each line read one int and three doubles    */
    /* corresponding to a vertex number and three coordinates of the vertex */
    this->ind_pos[0] = 0;
    for (int i=0; i < this->num_vertices; ++i){
        fscanf(fp, "%d %lf %lf %lf", &ivert,
                       &this->coordinates[3*i],
                       &this->coordinates[3*i+1],
                       &this->coordinates[3*i+2] );

        this->ind_pos[i+1] = this->ind_pos[i] + 3;
    }

    fclose(fp);
}

/******************************************************************************/

/* Makes a copy */
void meshgeom_copy(MeshGeometry* this, const MeshGeometry* other){

    int num_vertices = other->num_vertices;
    int num_coordinates = other->num_coordinates;

    this->num_vertices = num_vertices;
    this->num_coordinates = num_coordinates;

    this->ind_pos = realloc(this->ind_pos, (num_vertices+1)*sizeof(int));
    this->coordinates = realloc(this->coordinates,
            num_coordinates*sizeof(double));

    memcpy(this->ind_pos, other->ind_pos, (num_vertices+1)*sizeof(int));
    memcpy(this->coordinates, other->coordinates,
            num_coordinates*sizeof(double));
}

/******************************************************************************/

/* Deletes a MeshGeometry object */
void meshgeom_delete(MeshGeometry* this){

    free(this->ind_pos);
    free(this->coordinates);
    this->num_vertices = 0;
    this->num_coordinates = 0;
}

/******************************************************************************/

/* Clears a MeshGeometry object. */
void meshgeom_clear(MeshGeometry* this){

    free(this->ind_pos);
    free(this->coordinates);

    this->ind_pos = NULL;
    this->coordinates = NULL;

    this->num_vertices = 0;
    this->num_coordinates = 0;
}

/******************************************************************************/

/* Gets the number of vertices */
int meshgeom_get_num_vertices(const MeshGeometry* this){

    return this->num_vertices;
}

/******************************************************************************/

/* Gets the total number of coordinates */
int meshgeom_get_num_coordinates(const MeshGeometry* this){

    return this->num_coordinates;
}

/******************************************************************************/

/* Get a pointer to all coordinates */
double* meshgeom_get_coordinates(const MeshGeometry* this, int* num_coordinates){

    if (num_coordinates != NULL){
        *num_coordinates = this->num_coordinates;
    }
    return this->coordinates;
}

/******************************************************************************/

/* Gets a pointer to the coordinates of vertex ivert */
double* meshgeom_get_coords(const MeshGeometry* this, const int ivert, int* nc){

    if (nc != NULL){
        /* Three coordinates per vertex */
        *nc = 3;
    }
    return (this->coordinates + 3*ivert);
}

/******************************************************************************/

/* Gets a pointer to the coordinate i of vertex ivert; 0<=i<=2 */
double* meshgeom_get_coord(const MeshGeometry* this, const int ivert,
        const int i){

    return (this->coordinates + 3*ivert + i);
}

/******************************************************************************/

/* Sets all coordinates */
void meshgeom_set_coordinates(MeshGeometry* this, const double* coordinates){

    memcpy(this->coordinates, coordinates,
            this->num_coordinates*sizeof(double));
}

/******************************************************************************/

/* Sets the coordinates of vertex ivert */
void meshgeom_set_coords(MeshGeometry* this, const int ivert,
        const double* coords){

    memcpy(this->coordinates+3*ivert, coords, 3*sizeof(double));
}

/******************************************************************************/

/* Sets coordinate i of vertex ivert; 0<=i<=2 */
void meshgeom_set_coord(MeshGeometry* this, const int ivert, const int i,
        const double* coord){

    this->coordinates[this->ind_pos[ivert]+i] = *coord;
}

/******************************************************************************/

/* Writes to a text file */
void meshgeom_to_file(const MeshGeometry* this, const char* fn){

    FILE* fp = fopen(fn, "w");

    fprintf(fp, "%d  %d\n", this->num_vertices, this->num_coordinates);

    for (int i=0; i < this->num_vertices; ++i){
        fprintf(fp, "%d  %.17f  %.17f  %.17f\n", i,
                this->coordinates[3*i],
                this->coordinates[3*i+1],
                this->coordinates[3*i+2]);
    }

    fclose(fp);
}

/******************************************************************************/
