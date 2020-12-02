#include "ref_config.h"
#include "triangle_cell.h"
#include "utils_math.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/******************************************************************************/

RefConfig refconf_create()
{
    RefConfig rc;

    rc.num_cells = 0;
    rc.num_verts_per_cell = 0;
    rc.vertices = NULL;
    rc.normals = NULL;
    rc.A = NULL;
    rc.B = NULL;

    return rc;
}

/******************************************************************************/

void refconf_delete(RefConfig* this)
{
    free(this->vertices);
    free(this->normals);
    free(this->A);
    free(this->B);
}

/******************************************************************************/

void refconf_init(RefConfig* this, Mesh* msh)
{
    /* TODO: Fix the condition for cell_type */
    if (strcmp(msh->cell_type, "tri")==0){
        this->num_verts_per_cell = 3;
    }
    else{
        exit(EXIT_FAILURE);
    }

    this->num_cells = mesh_get_num_entities(msh, 2);

    /* x_0 = y_0 = y_1 = 0, hence 3 coordinates per triangle */
    this->vertices = malloc(3*this->num_cells*sizeof(double));
    this->normals = malloc(3*this->num_cells*sizeof(double));
    this->A = malloc(3*this->num_cells*sizeof(double));
    this->B = malloc(3*this->num_cells*sizeof(double));

    /* Calculating the vertex coordinates in local basis for each triangular cell */
    int nv;
    int nc;
    int* icell_verts;
    double* coords; 
    double L;
    double icell_coords[9]; /* number of local coordinates */
    double tmp[9];
    double x[3], y[3], a[3], b[3], c[3];

    for(int icell=0; icell < this->num_cells; ++icell){
        icell_verts = mesh_get_conns(msh, 2, 0, icell, &nv);

        for(int i=0; i< nv; ++i){
            coords = mesh_get_coords(msh, icell_verts[i], &nc);
            memcpy((icell_coords+3*i), coords, nc*sizeof(double));
        }

        tri_to_local(icell_coords, tmp, NULL);

        //printf("tmp:ref_config:icell=%d\n", icell);
        //for (int i=0; i < 3; ++i){
        //    for (int j=0; j < 3; ++j){
        //        printf("%f  ", tmp[3*i+j]);
        //    }
        //    printf("\n");
        //}

        this->vertices[3*icell] = tmp[3];
        this->vertices[3*icell+1] = tmp[6];
        this->vertices[3*icell+2] = tmp[7];

        tri_normal(icell_coords, this->normals+3*icell);

        for (int i=0; i < 3; ++i){
            x[i] = tmp[3*i+0];
            y[i] = tmp[3*i+1];
        }

        a[0] = y[1] - y[2];
        a[1] = y[2] - y[0];
        a[2] = y[0] - y[1];

        b[0] = x[2] - x[1];
        b[1] = x[0] - x[2];
        b[2] = x[1] - x[0];

        cross(x, y, c);

        for (int i=0; i < 3; ++i){
            L = a[i]*x[i] + b[i]*y[i] + c[i];
            this->A[3*icell+i] = a[i]/L;
            this->B[3*icell+i] = b[i]/L;
        }
    }
}

/******************************************************************************/

/* Returns a pointer to an array containing y1, x2, and y2 */
double* refconf_get_coords(RefConfig* this, int icell, int* nc)
{
    if (nc != NULL){
        *nc = 3;
    }
    return (this->vertices+3*icell);
}

/******************************************************************************/

/* Returns a pointer to vector A for cell icell */
double* refconf_get_A(RefConfig* this, int icell, int* nc)
{
    if (nc != NULL){
        *nc = 3;
    }
    return (this->A+3*icell);
}

/******************************************************************************/

/* Returns a pointer to vector B for cell icell */
double* refconf_get_B(RefConfig* this, int icell, int* nc)
{
    if (nc != NULL){
        *nc = 3;
    }
    return (this->B+3*icell);
}

/******************************************************************************/
