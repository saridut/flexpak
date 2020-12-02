#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mesh.h"
#include "boundary_mesh.h"
#include "ref_config.h"

struct Assembler{
    Mesh* msh;
    BoundaryMesh* bmsh;
    RefConfig* refconf;
    double* cell_center;
    double* cell_area;
    double* cell_angle;
    double* cell_angle_cot;
    double* cell_normal;
    double* edge_len_sq;
    double* edge_center;
    double* amixed;
    double* vert_normal;
    double* mean_curv;
    double* gauss_curv;
    double* sl_op;
    double* sgmc;
    double* slmc;
    int* is_obtuse;
    double* traction;
    double E;
    double kb;

};

typedef struct Assembler Assembler;


Assembler asm_create();

/* Allocate memory and setup buffers for force assembly. If bending */
/* effects are not required pass a negative value for kb. If there  */
/* is no boundary, pass a null pointer as bmsh.                     */
void asm_init(Assembler* this, Mesh* msh, BoundaryMesh* bmsh,
        RefConfig* refconf, double E, double kb);


/* Calculates the single layer density q at all nodes */
void asm_calc_sldq(Assembler* this);


/* Calculates the traction at all nodes */
void asm_calc_traction(Assembler* this);


/* Returns a pointer to the traction vector at all nodes */
double* asm_get_traction(const Assembler* this, int* n);


/* Returns a pointer to the traction vector at vert ivert */
double* asm_get_traction_at_vert(const Assembler* this, const int ivert);


/* Returns a pointer to the single layer density q at all nodes */
double* asm_get_sldq(const Assembler* this, int* n);


void asm_delete(Assembler* this);


void asm_calc_traction_tension(Assembler* this);


void asm_calc_traction_bending(Assembler* this);

void asm_set_bc(Assembler* this);

void asm_calc_edge_len_sq(Assembler* this);

void asm_calc_angles(Assembler* this);

void asm_calc_edge_centers(Assembler* this);

void asm_calc_cell_centers(Assembler* this);

void asm_calc_sgmc(Assembler* this);


#ifdef __cplusplus
}
#endif

#endif /* ASSEMBLER_H */
