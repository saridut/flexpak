#ifndef REF_CONFIG_H
#define REF_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mesh.h"

struct RefConfig{
    int num_cells;
    int num_verts_per_cell;
    double* vertices;
    double* normals;
    double* A;
    double* B;
};

typedef struct RefConfig RefConfig;


RefConfig refconf_create();


void refconf_delete(RefConfig* this);


void refconf_init(RefConfig* this, Mesh* msh);


/* Returns a pointer to an array containing y1, x2, and y2 */
double* refconf_get_coords(RefConfig* this, int icell, int* nc);


/* Returns a pointer to vector A for cell icell */
double* refconf_get_A(RefConfig* this, int icell, int* nc);


/* Returns a pointer to vector B for cell icell */
double* refconf_get_B(RefConfig* this, int icell, int* nc);

#ifdef __cplusplus
}
#endif

#endif /* REF_CONFIG_H */
