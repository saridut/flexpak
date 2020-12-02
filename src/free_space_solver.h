#ifndef FREE_SPACE_SOLVER_H
#define FREE_SPACE_SOLVER_H

#include "mesh.h"

struct FreeSpaceSolver{
    Mesh* msh;
    double* velocity;
    /*! Single layer density (denoted q in Amit's paper) */
    double* sldq; 
    double gam_dot;
    /*! Number of quadrature points */
    int num_quad_points[2];
};

typedef struct FreeSpaceSolver FreeSpaceSolver;


FreeSpaceSolver fs_solver_create();


void fs_solver_add_mesh(FreeSpaceSolver* this, Mesh* msh);


void fs_solver_set_shear_rate(FreeSpaceSolver* this, double gam_dot);

/*! Same in along both coordinate directions. Good to keep \p n1 = \p n2. */
void fs_solver_set_num_quad_points(FreeSpaceSolver* this, int n1, int n2);


void fs_solver_set_sldq(FreeSpaceSolver* this, double* q);


void fs_solver_init(FreeSpaceSolver* this);


void fs_solver_delete(FreeSpaceSolver* this);


void fs_solver_solve(FreeSpaceSolver* this);


void fs_solver_calc_sli(FreeSpaceSolver* this);


double* fs_solver_get_velocity(FreeSpaceSolver* this);


void fs_solver_add_ambient_velocity(FreeSpaceSolver* this);


double fs_solver_get_area_element(double x[3], double y[3]);


void fs_solver_calc_green_func(const double r[3], double* green_func);


void fs_solver_calc_sli_nosing(FreeSpaceSolver* this, const double x0[3],
        const double x1[3], const double x2[3], const double x3[3],
        const double q1[3], const double q2[3], const double q3[3],
        double u[3]);


void fs_solver_calc_sli_sing(FreeSpaceSolver* this, const double x0[3],
        const double x1[3], const double x2[3], const double x3[3],
        const double q1[3], const double q2[3], const double q3[3],
        double u[3]);


#endif /* FREE_SPACE_SOLVER_H */

