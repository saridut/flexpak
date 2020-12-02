#ifndef GGEM_H
#define GGEM_H

#include "mesh.h"
#include "local_cells.hpp"
#include <complex.h>
#include <fftw3.h>

#define Lx 30.0
#define Ly 30.0
#define Lz 30.0

#define nx 128
#define ny 128
#define nz 129

struct GgemSolver{
    Mesh* msh;
    LocalCells lc;
    double alpha;
    double rcut;
    double* velocity;
    double* sldq;
    double* u_gl;
    double* v_gl;
    double* w_gl;
    double* p_gl;
    /* Used as global variables in Amit's code */

    // chebyshev points and normalization coefficients
    // following are pointers to array of size nz
    double xz[nz];
    double xxi[nz];
    double cz[nz];

    //Solutions to the two homogeneous problems
    double complex u1H[nx][ny/2+1][nz];
    double complex v1H[nx][ny/2+1][nz];
    double complex w1H[nx][ny/2+1][nz];
    double complex p1H[nx][ny/2+1][nz];
    double complex dudz1H[nx][ny/2+1][nz];
    double complex dvdz1H[nx][ny/2+1][nz];
    double complex dwdz1H[nx][ny/2+1][nz];

    double complex* u2H[nx][ny/2+1][nz];
    double complex* v2H[nx][ny/2+1][nz];
    double complex* w2H[nx][ny/2+1][nz];
    double complex* p2H[nx][ny/2+1][nz];
    double complex* dudz2H[nx][ny/2+1][nz];
    double complex* dvdz2H[nx][ny/2+1][nz];
    double complex* dwdz2H[nx][ny/2+1][nz];
};

typedef struct GgemSolver GgemSolver;


GgemSolver ggem_solver_create();


void ggem_add_mesh(GgemSolver* this, Mesh* msh);


void ggem_set_sldq(GgemSolver* this, double* q)


void ggem_init(GgemSolver* gs);


void ggem_solve_init(GgemSolver* gs);


void ggem_solve(GgemSolver* gs);


double* ggem_get_velocity(GgemSolver* gs);


void ggem_solver_delete(GgemSolver* gs);


void ggem_calc_sli(GgemSolver* gs, double* q);


void ggem_calc_sli_local(const double x0[3], const double x1[3],
        const double x2[3], const double x3[3], const double f1[3], 
        const double f2[3], const double f3[3], const double alpha, double u[3]);


void ggem_calc_sli_local_polar(const double x0[3], const double x1[3],
        const double x2[3], const double x3[3], const double f1[3], 
        const double f2[3], const double f3[3], const double alpha, double u[3]);


void ggem_sldq_to_glmesh(double xb[beads][3], double f[beads][3],
        double gaussx[nx][ny][nz], double gaussy[nx][ny][nz],
        double gaussz[nx][ny][nz], double gaussz_p[nx][ny][nz],
        double weight[beads], double alpha, double xz[nz]);


/* Interpolate the global velocity at the field point                  */
/* Uses quartic interpolation between the nodes of the containing cell */
/* currently assumes that the points are well inside the walls         */
void ggem_velocity_from_glmesh(double xb[beads][3], double ub[beads][3],
        double uxf[nx][ny][nz], double uyf[nx][ny][nz],
        double uzf[nx][ny][nz], double xz[nz]);


void solve_quasi_tridiagonal(double nu, double b,double complex fhat[nz],
            double complex cheb_bc[2], double complex uhat[nz], double cz[nz]);


void compute_coef_derivative(double complex u[nz],double complex uprime[nz],
            double cz[nz]);

void ggem_global_velocity_homogeneous(
        double complex rhox[nx][ny/2+1][nz], 
        double complex rhoy[nx][ny/2+1][nz],
        double complex rhoz[nx][ny/2+1][nz],
        double complex rhoz_p[nx][ny/2+1][nz],
	    double complex dudz[nx][ny/2+1][nz],
        double complex dvdz[nx][ny/2+1][nz],
	    double complex ptr[nx][ny/2+1][nz],
        double cz[nz],
        int problem_no);


void ggem_global_velocity_inhomogeneous(
        double gaussx[nx][ny][nz], double gaussy[nx][ny][nz],
        double gaussz[nx][ny][nz], double gaussz_p[nx][ny][nz],
        double ubc1[nx][ny][3], double ubc2[nx][ny][3],
        double uxf[nx][ny][nz], double uyf[nx][ny][nz], double uzf[nx][ny][nz],
        double pf[nx][ny][nz], 
        double dudx[nx][ny][nz], double dudy[nx][ny][nz], double dudz[nx][ny][nz],
        double dvdx[nx][ny][nz], double dvdy[nx][ny][nz], double dvdz[nx][ny][nz],
        double dwdx[nx][ny][nz], double dwdy[nx][ny][nz], double dwdz[nx][ny][nz],
        double xxi[nz], double cz[nz],
        double complex rhox1[nx][ny/2+1][nz], double complex rhoy1[nx][ny/2+1][nz],
        double complex rhoz1[nx][ny/2+1][nz], double complex rhoz_p1[nx][ny/2+1][nz],
        double complex dudz1[nx][ny/2+1][nz], double complex dvdz1[nx][ny/2+1][nz], 
        double complex ptr1[nx][ny/2+1][nz], 							 
        double complex rhox2[nx][ny/2+1][nz], double complex rhoy2[nx][ny/2+1][nz],
        double complex rhoz2[nx][ny/2+1][nz], double complex rhoz_p2[nx][ny/2+1][nz],
        double complex dudz2[nx][ny/2+1][nz], double complex dvdz2[nx][ny/2+1][nz],
        double complex ptr2[nx][ny/2+1][nz]);







#endif /* GGEM_H */
