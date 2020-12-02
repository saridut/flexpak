#ifndef MESH_GEOMETRY_H
#define MESH_GEOMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

struct MeshGeometry{
    int num_vertices;
    int num_coordinates;
    double* coordinates;
    int* ind_pos;
};

typedef struct MeshGeometry MeshGeometry;

/* Creates an empty MeshGeometry object */
MeshGeometry meshgeom_create();

/* Initializes from an array of coordinates.*/
void meshgeom_from_array(MeshGeometry* this, const int num_vertices,
        const double* coordinates);

/* Initializes from file */
void meshgeom_from_file(MeshGeometry* this, const char* fn);

/* Makes a copy */
void meshgeom_copy(MeshGeometry* this, const MeshGeometry* other);

/* Deletes a MeshGeometry object */
void meshgeom_delete(MeshGeometry* this);

/* Clears a MeshGeometry object. */
void meshgeom_clear(MeshGeometry* this);

/* Gets the number of vertices */
int meshgeom_get_num_vertices(const MeshGeometry* this);

/* Gets the total number of coordinates */
int meshgeom_get_num_coordinates(const MeshGeometry* this);

/* Gets a pointer to all coordinates */
double* meshgeom_get_coordinates(const MeshGeometry* this,
        int* num_coordinates);

/* Gets a pointer to the coordinates of vertex ivert */
double* meshgeom_get_coords(const MeshGeometry* this, const int ivert,
        int* nc);

/* Gets a pointer to the coordinate i of vertex ivert; 0<=i<=2 */
double* meshgeom_get_coord(const MeshGeometry* this, const int ivert,
        const int i);

/* Sets all coordinates */
void meshgeom_set_coordinates(MeshGeometry* this, const double* coordinates);

/* Sets the coordinates of vertex ivert */
void meshgeom_set_coords(MeshGeometry* this, const int ivert,
        const double* coords);

/* Sets coordinate i of vertex ivert; 0<=i<=2 */
void meshgeom_set_coord(MeshGeometry* this, const int ivert, const int i,
        const double* coord);

/* Writes to a text file */
void meshgeom_to_file(const MeshGeometry* this, const char* fn);

#ifdef __cplusplus
}
#endif

#endif //MESH_GEOMETRY_H
