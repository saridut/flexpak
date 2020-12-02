#ifndef TRIANGLE_CELL_H
#define TRIANGLE_CELL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Calculates the area of a triangle */
double tri_area(const double* coords);


/* Calculates the normal to a triangle */
void tri_normal(const double* coords, double* normal);


/* Find the circumcenter of a triangle */
void tri_circumcenter(double* coords, double* cc, double* cr,
        double* xi, double* eta);


/* Calculates the projection operator `op' for a given normal. */
/* op = I - pp, where p is the normal vector.                  */
void op_proj(const double* normal, double* op);


/* Calculates coordinates of vertices in local basis */
void tri_to_local(const double* g_verts, double* l_verts, double* dcm);


void tri_get_local_basis(const double* verts, double* basis);


void tri_get_displacement(const double* verts_A, const double* verts_B,
        double* u, double* v);

#ifdef __cplusplus
}
#endif

#endif /* TRIANGLE_CELL_H */
