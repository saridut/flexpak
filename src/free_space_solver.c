#include "free_space_solver.h"
#include "utils_math.h"
#include "gauss_legendre.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/******************************************************************************/

FreeSpaceSolver fs_solver_create()
{
    FreeSpaceSolver fss;

    fss.msh = NULL;
    fss.velocity = NULL;
    fss.sldq = NULL;
    fss.gam_dot = 0.0;
    fss.num_quad_points[0] = 0;
    fss.num_quad_points[1] = 0;

    return fss;
}

/******************************************************************************/

void fs_solver_add_mesh(FreeSpaceSolver* this, Mesh* msh)
{
    this->msh = msh;
}

/******************************************************************************/

void fs_solver_set_shear_rate(FreeSpaceSolver* this, double gam_dot)
{
    this->gam_dot = gam_dot;
}

/******************************************************************************/

void fs_solver_set_num_quad_points(FreeSpaceSolver* this, int n1, int n2)
{
    this->num_quad_points[0] = n1;
    this->num_quad_points[1] = n2;
}

/******************************************************************************/

void fs_solver_set_sldq(FreeSpaceSolver* this, double* q)
{
    this->sldq = q;
}

/******************************************************************************/

void fs_solver_init(FreeSpaceSolver* this)
{
    int num_coordinates = mesh_get_num_coordinates(this->msh);
    this->velocity = malloc(num_coordinates*sizeof(double));
}

/******************************************************************************/

void fs_solver_delete(FreeSpaceSolver* this)
{
    this->msh = NULL;
    this->sldq = NULL;
    free(this->velocity);
}

/******************************************************************************/

void fs_solver_solve(FreeSpaceSolver* this)
{
    int num_coordinates = mesh_get_num_coordinates(this->msh);

    fs_solver_calc_sli(this);

    fs_solver_add_ambient_velocity(this);
}

/******************************************************************************/

double* fs_solver_get_velocity(FreeSpaceSolver* this)
{
    return this->velocity;
}

/******************************************************************************/

void fs_solver_calc_sli(FreeSpaceSolver* this)
{
    double u[3];
    double  *x0, *x1, *x2;
    double  *q0, *q1, *q2;
    double* coords_ivert;
    int* conns;
    int nc;
    int conns_buf[3]; //For tri cells, size 3

    int num_verts = mesh_get_num_vertices(this->msh);
    int num_cells = mesh_get_num_entities(this->msh, 2);

    zero_out(3*num_verts, this->velocity);

    /* Contribution from local force density */
    for(int ivert=0; ivert < num_verts; ++ivert) {
        coords_ivert = mesh_get_coords(this->msh, ivert, NULL);
        for (int icell=0; icell < num_cells; ++icell){
            conns = mesh_get_conns(this->msh, 2, 0, icell, &nc);
            for (int i=0; i < nc; ++i){
                conns_buf[i] = conns[i];
            }
            /* Permute vertices */
            if (ivert==conns_buf[1]){
                conns_buf[1] = conns_buf[2];
                conns_buf[2] = conns_buf[0];
                conns_buf[0] = ivert;
            }
            else if (ivert==conns[2]){
                conns_buf[2] = conns_buf[1];
                conns_buf[1] = conns_buf[0];
                conns_buf[0] = ivert;
            }
            x0 = mesh_get_coords(this->msh, conns_buf[0], NULL);
            x1 = mesh_get_coords(this->msh, conns_buf[1], NULL);
            x2 = mesh_get_coords(this->msh, conns_buf[2], NULL);

            q0 = &this->sldq[3*conns_buf[0]];
            q1 = &this->sldq[3*conns_buf[1]];
            q2 = &this->sldq[3*conns_buf[2]];

            //for(int j=0; j < 3; ++j){
            //    printf("q0: %f  %f  %f\n", q0[0], q0[1], q0[2]);
            //    printf("q1: %f  %f  %f\n", q1[0], q1[1], q1[2]);
            //    printf("q2: %f  %f  %f\n", q2[0], q2[1], q2[2]);
            //}
            if (mesh_is_incident(this->msh, 0, ivert, 2, icell)){
                fs_solver_calc_sli_sing(this, coords_ivert, x0, x1, x2, q0, q1, q2, u);
            }
            else {
                fs_solver_calc_sli_nosing(this, coords_ivert, x0, x1, x2, q0, q1, q2, u);
            }
            for (int j=0; j < 3; ++j){
                this->velocity[3*ivert+j] += u[j];
            }
        }
    }
}

/******************************************************************************/

//Only sinple shear flow implemented: vx = gam_dot*y
void fs_solver_add_ambient_velocity(FreeSpaceSolver* this)
{
    int num_verts = mesh_get_num_vertices(this->msh);
    double* coordinates = mesh_get_coordinates(this->msh, NULL);

    for (int ivert=0; ivert < num_verts; ++ivert){
        this->velocity[3*ivert] += this->gam_dot*coordinates[3*ivert+1];
    }
}

/******************************************************************************/

/*! Calculates the area element (required for evaluating a surface integral) */
double fs_solver_get_area_element(double x[3], double y[3])
{
    double z[3];

    cross(x, y, z);
    return norm(3, z);
}

/******************************************************************************/

/*! Calculates the Green function for an unbounded domain. Note that the
 * prefactor of 1/8\pi is not included. */
void fs_solver_calc_green_func(const double r[3], double* green_func)
{
    const int ncols = 3;
    double rg  = norm(3, r);
    double rg_inv = 1.0/rg;
    double rg_inv_sq = rg_inv*rg_inv;
    double rg_inv_cu = rg_inv*rg_inv_sq;
    
    /* find the local green's function */
    green_func[ncols*0+0] = rg_inv+rg_inv_cu*r[0]*r[0];
    green_func[ncols*1+1] = rg_inv+rg_inv_cu*r[1]*r[1];
    green_func[ncols*2+2] = rg_inv+rg_inv_cu*r[2]*r[2];
    
    green_func[ncols*0+1] = rg_inv_cu*r[0]*r[1];
    green_func[ncols*0+2] = rg_inv_cu*r[0]*r[2];
    green_func[ncols*1+2] = rg_inv_cu*r[1]*r[2];
    
    green_func[ncols*1+0] = green_func[3*0+1];
    green_func[ncols*2+0] = green_func[3*0+2];
    green_func[ncols*2+1] = green_func[3*1+2];
}

/******************************************************************************/

/*!
 *  @brief Evaluates the single layer integral over a triangle
 *
 *  This is the non-sigular version.
 *
 *  @param x0 Field point
 *  @param x1, x2, x3  Vertices of a triangle
 *  @param q1, q2, q3  Value of vector q at the triangle vertices
 *  @param u [out] Value of the integral at \p x0
 *
 *  @todo change the integration rule from tensor product basis
 *
 *  @warning Really no idea what is going on here
 */
void fs_solver_calc_sli_nosing(FreeSpaceSolver* this, const double x0[3],
        const double x1[3], const double x2[3], const double x3[3],
        const double q1[3], const double q2[3], const double q3[3],
        double u[3])
{
    int nphi = this->num_quad_points[0];
    int nr = this->num_quad_points[1];

    double* zphi = gleg_get_nodes(nphi);
    double* wphi = gleg_get_weights(nphi);
    double* zr = gleg_get_nodes(nr);
    double* wr = gleg_get_weights(nphi);

    /* find the area element omega; constant for linear element */
    double dx_dxi[3];
    double dx_deta[3];

    for(int i=0; i < 3; ++i) {
        dx_dxi[i] = x3[i] - x1[i];
        dx_deta[i] = x2[i] - x1[i];
    }
    double omega = fs_solver_get_area_element(dx_dxi, dx_deta);

    zero_out(3, u);

    double l1, l2, l3;
    double xi, eta;
    double q[3], x[3], rhat[3];
    double green_func[9];
    double w;

    for(int iphi=0; iphi < nphi; ++iphi) {
        eta = 0.5*(1.0+zphi[iphi]);

        for(int ir=0; ir < nr; ++ir) {
            xi = 0.5*(1.0+zr[ir])*(1.0-eta);

            /* find the xi and eta coordinates */
            l1 = 1.0-xi-eta;
            l2 = eta;
            l3 = xi;

            /* coefficient of green's function at the current point */
            for(int i=0; i<3; ++i) {
                q[i] = l1*q1[i] + l2*q2[i] + l3*q3[i];
            }

            /* position of current point */
            for(int i=0; i < 3; ++i) {
                x[i] = l1*x1[i] + l2*x2[i] + l3*x3[i];
            }

            /* separation vector: field - pole */
            for(int i=0; i < 3; ++i) {
                rhat[i] = x0[i] - x[i];
            }

            /* global coordinates distance rg */
            fs_solver_calc_green_func(rhat, green_func);

            /* weight */
            w  = wphi[iphi]*wr[ir]*(1/2.0)*(1-eta)/2.0*omega;

            /* find contribution to velocity */
            for(int i=0; i < 3; ++i) {
                for(int j=0; j < 3; ++j) {
                    u[i] += q[j]*green_func[3*i+j]*w;
                }
            }

        }
    }
}

/******************************************************************************/

/*!
 *  @brief Evaluates the single layer integral over a triangle
 *
 *  This is the sigular version. The singularity is handled by transforming to
 *  polar coordinates.
 *
 *  @param x0 Field point
 *  @param x1, x2, x3  Vertices of a triangle
 *  @param q1, q2, q3  Value of vector q at the triangle vertices
 *  @param u [out] Value of the integral at \p x0
 *
 *  @warning Really no idea what is going on here
 */
void fs_solver_calc_sli_sing(FreeSpaceSolver* this, const double x0[3],
        const double x1[3], const double x2[3], const double x3[3],
        const double q1[3], const double q2[3], const double q3[3],
        double u[3])
{
    int nphi = this->num_quad_points[0];
    int nr = this->num_quad_points[1];

    double* zphi = gleg_get_nodes(nphi);
    double* wphi = gleg_get_weights(nphi);
    double* zr = gleg_get_nodes(nr);
    double* wr = gleg_get_weights(nphi);

    /* find the area element omega */
    double dx_dxi[3];
    double dx_deta[3];

    for(int i=0; i < 3; ++i) {
        dx_dxi[i]  = x3[i] - x1[i];
        dx_deta[i] = x2[i] - x1[i];
    }
    /* area element is constant for linear element */
    double omega = fs_solver_get_area_element(dx_dxi, dx_deta); 

    zero_out(3, u);
    
    double phi1 = 0.0;
    double phi2 = M_PI_2;
    double r1, r2;
    double phi, r;
    double xi, eta;
    double l1, l2, l3;
    double q[3], x[3], rhat[3];
    double green_func[9];
    double w;

    for(int iphi=0; iphi < nphi; ++iphi) {
        phi = phi1 + 0.5*(1.0+zphi[iphi])*(phi2-phi1);
        r1 = 0.0;
        r2 = sin(M_PI/4.0)/sin(3*M_PI/4.0-phi);

        for(int ir = 0; ir < nr; ++ir) {
            r = r1 + 0.5*(1.0+zr[ir])*(r2-r1);

            /* find the xi and eta coordinates */
            xi = r*sin(phi);
            eta = r*cos(phi);

            /* find the xi and eta coordinates */
            l1 = 1.0-xi-eta;
            l2 = eta;
            l3 = xi;

            /* coefficient of green's function at the current point */
            for(int i=0; i<3; ++i) {
                q[i] = l1*q1[i] + l2*q2[i] + l3*q3[i];
            }

            /* position of current point */
            for(int i=0; i<3; ++i) {
                x[i] = l1*x1[i]+ l2*x2[i] + l3*x3[i];
            }

            /* separation vector: field - pole */
            for(int i=0; i<3; ++i) {
                rhat[i] = x0[i] - x[i];
            }

            fs_solver_calc_green_func(rhat, green_func);

            /* weight */
            w  = wphi[iphi]*wr[ir]*(M_PI/2.0/2.0)*(r2-r1)/2.0*omega;

            /* find contribution to velocity */
            for(int i=0; i<3; ++i) {
                for(int j=0; j<3; ++j) {
                    u[i] += q[j]*green_func[3*i+j]*w*r;
                }
            }
        }
    }
}

/******************************************************************************/

