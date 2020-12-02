#include "ggem.h"
#include <stddef.h>
#include <math.h>

/******************************************************************************/

GgemSolver ggem_solver_create()
{
    GgemSolver gs;

    gs.msh = NULL;
    gs.velocity = NULL;
    gs.sldq = NULL;
    gs.u_gl = NULL;
    gs.v_gl = NULL;
    gs.w_gl = NULL;
    gs.p_gl = NULL;

    return gs;
}

/******************************************************************************/

void ggem_add_mesh(GgemSolver* this, Mesh* msh)
{
    this->msh = msh;
}

/******************************************************************************/

void ggem_set_sldq(GgemSolver* this, double* q)
{
    this->sldq = q;
}

/******************************************************************************/

void ggem_init(GgemSolver* this)
{
    this->rcut = 4.0/this->alpha;
    int num_verts = mesh_get_num_vertices(this->msh);
    this->velocity = malloc(3*num_verts);

    //Solutions to the homogeneous problems are zeroed out in the
    //function ggem_global_velocity_homogeneous.

    //Setting Chebyshev nodes -- xxi in [-1, 1] & xz in [0, Lz]
    //Chebyshev coefficients -- cz
    for (int i=0; i<nz; ++i) {
		this->xxi[i] = -cos(M_PI*i/(nz-1));
		this->xz[i] = Lz/2.0*(1.0+xxi[i]);
        this->cz[i] = 1.0;
    }
    this->cz[0] = 2.0;
    this->cz[nz-1] = 2.0;

    /* Set integration points and weights */
    set_pos_weights_singular()
    set_pos_weights_nonsingular()


}

/******************************************************************************/

void ggem_solver_delete(GgemSolver* this)
{
    this->msh = NULL;
    this->sldq = NULL;
    free(this->velocity);
}

/******************************************************************************/

void ggem_solve_init(GgemSolver* this)
{
    /* solve the homogeneous global solutions */
	ggem_global_velocity_homogeneous(this->u1H, this->v1H, this->w1H,
            this->dwdz1H, this->dudz1H, this->dvdz1H,
            this->p1H, this->cz, 1);

	ggem_global_velocity_homogeneous(this->u2H, this->v2H, this->w2H,
            this->dwdz2H, this->dudz2H, this->dvdz2H,
            this->p2H, this->cz, 2);
}

/******************************************************************************/

void ggem_solve(GgemSolver* this)
{
    ggem_calc_sli(this);
    ggem_add_ambient_velocity(this);
}

/******************************************************************************/

double* ggem_get_velocity(GgemSolver* this)
{
    return this->velocity;
}

/******************************************************************************/

void ggem_calc_sli(GgemSolver* this)
{
    int i,j,k;
    int ibead,itriangle,jtriangle;
    double  xi[3];
    double r1;
    int i1,i2,i3;
    double pi;
    double u[3];
    int member;

    double weights[mesh_get_num_vertices(this->msh)];

    double uxf[nx][ny][nz];
    double uyf[nx][ny][nz];
    double uzf[nx][ny][nz];
    double pf[nx][ny][nz];

    double ubc1[nx][ny][3];
    double ubc2[nx][ny][3];

    double gaussx[nx][ny][nz];
    double gaussy[nx][ny][nz];
    double gaussz[nx][ny][nz];
    double gaussz_p[nx][ny][nz];

    int ix,iy,iz;
    double dx,dy;
    int ipart,jpart;
    double fn,f1[3],f2[3],f3[3];
    double* coords_ivert;
    int* conns;
    int* loc_cells;
    int nc;
    int nloc;

    /* Mesh size */
    dx = Lx/nx;
    dy = Ly/ny;

    /* Zeroing out velocity */
    int num_verts = mesh_get_num_vertices(this->msh);
    for (int ivert=0; ivert < num_verts; ++ivert) {
        for (int j=0; j < 3; ++j) {
            this->velocity[3*ivert+j] = 0.0;
        }
    }

    /* Contribution from local force density */
    for(int ivert=0; ivert < num_verts; ++ivert) {
        coords_ivert = mesh_get_coords(this->msh, ivert, NULL);
        loc_cells = (this->lc).get_cells(ivert, &nloc);
        for (int icell=0; icell < nloc; ++icell){
            conns = mesh_get_conns(this->msh, 2, 0, icell, &nc);
            if (ivert==conns[1]){
                conns[1] = conns[2];
                conns[2] = conns[0];
                conns[0] = ivert;
            }
            else if (ivert==conns[2]){
                conns[2] = conns[1];
                conns[1] = conns[0];
                conns[0] = ivert;
            }
            x0 = mesh_get_coords(this->msh, conns[0], NULL);
            x1 = mesh_get_coords(this->msh, conns[1], NULL);
            x2 = mesh_get_coords(this->msh, conns[2], NULL);

            if (mesh_is_incident(this->msh, 0, ivert, 2, icell)){
                ggem_calc_sli_local_polar(coords_ivert, x0, x1, x2,
                    &sldq[3*conns[0]], &sldq[3*conns[1]], &sldq[3*conns[2]],
                    this->alpha, u);
            }
            else {
                ggem_calc_sli_local(coords_ivert, x0, x1, x2,
                    &sldq[3*conns[0]], &sldq[3*conns[1]], &sldq[3*conns[2]],
                    this->alpha, u);
            }

            for (int j=0; j < 3; ++j){
                this->velocity[3*ivert+j] += u[j];
            }

        }
    }

    /*---------------- Global solution calculation -----------------*/

    /*----- Boundary Conditions ---- */
    /* Zeroing out */
    for (int i=0; i < nx; ++i) {
        for (int j=0; j < ny; ++j) {
            for (int k=0; k < 3; ++k) {
                ubc1[i][j][k] = 0.0;
                ubc2[i][j][k] = 0.0;
            }
        }
    }

    /* set ubc = - u_l for poiseuille flow (??) */
    /* Top plate */
    for (int ix=0; ix<nx; ++ix) {
        for (int iy=0; iy<ny; ++iy) {
            xi[0] = ix*dx;
            xi[1] = iy*dy;
            xi[2] = Lz;

            for (int ivert=0; ivert < num_verts; ++ivert) {
                coords_ivert = mesh_get_coords(this->msh, ivert, NULL);
                loc_cells = (this->lc).get_cells(ivert, &nloc);
                for (int icell=0; icell < nloc; ++icell){
                    conns = mesh_get_conns(this->msh, 2, 0, icell, &nc);
                    if (ivert==conns[1]){
                        conns[1] = conns[2];
                        conns[2] = conns[0];
                        conns[0] = ivert;
                    }
                    else if (ivert==conns[2]){
                        conns[2] = conns[1];
                        conns[1] = conns[0];
                        conns[0] = ivert;
                    }
                    x0 = mesh_get_coords(this->msh, conns[0], NULL);
                    x1 = mesh_get_coords(this->msh, conns[1], NULL);
                    x2 = mesh_get_coords(this->msh, conns[2], NULL);

            for(itriangle=0; itriangle<triangle_count_t[ix][iy]; itriangle++) {
                jpart = triangle_list_t[ix][iy][itriangle][0];
                jtriangle = triangle_list_t[ix][iy][itriangle][1];

                i1 = triangle[jtriangle][0];
                i2 = triangle[jtriangle][1];
                i3 = triangle[jtriangle][2];
                member=0; // performing non-singular integral as of now

                ggem_calc_sli_local(xi, x0, x1, x2, &sldq[3*conns[0]],
                        &sldq[3*conns[1]], &sldq[3*conns[2]], this->alpha, u);

                integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],
                               fb[jpart][i2],fb[jpart][i3],u,alpha);


                for(int i=0; i<3; ++i) {
                    ubc2[ix][iy][i] -= u[i]/8/M_PI;
                }
            }
        }
    }

    /* bottom plate */
    for(ix=0; ix<nx; ix++) {
        for(iy=0; iy<ny; iy++) {
            xi[0] = ix*dx;
            xi[1] = iy*dy;
            xi[2] = 0.0;

            for(itriangle=0; itriangle<triangle_count_b[ix][iy]; itriangle++) {
                jpart = triangle_list_b[ix][iy][itriangle][0];
                jtriangle = triangle_list_b[ix][iy][itriangle][1];

                i1 = triangle[jtriangle][0];
                i2 = triangle[jtriangle][1];
                i3 = triangle[jtriangle][2];
                member=0; // performing non-singular integral as of now

                if (member == 1) {
                    integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],
                               fb[jpart][i2],fb[jpart][i3],u,alpha);
                }
                else {
                    integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],
                               fb[jpart][i2],fb[jpart][i3],u,alpha);
                }


                for(i=0; i<3; i++) {
                    ubc1[ix][iy][i] -= u[i]/8/pi;
                }

            }
        }
    }

    /* Find weights for the global solution */
    /* Initialize weights for the global solution */
    for(ipart=0; ipart<npart; ipart++) {
        for(i=0; i<beads; i++) {
            weights[ipart][i]=0.0;
        }
    }
    for(ipart=0; ipart<npart; ipart++) {
        for(itriangle=0; itriangle<ntriangles; itriangle++) {
            i1 = triangle[itriangle][0];
            i2 = triangle[itriangle][1];
            i3 = triangle[itriangle][2];

            /* Store weights for the global calculation */
            weights[ipart][i1] += area[ipart][itriangle]/3.0;
            weights[ipart][i2] += area[ipart][itriangle]/3.0;
            weights[ipart][i3] += area[ipart][itriangle]/3.0;

        }
    }

    /* Distribute density to mesh */
    /* Initialize mesh force density */
    for(ix=0; ix<nx; ix++) {
        for(iy=0; iy<ny; iy++) {
            for(iz=0; iz<nz; iz++) {
                gaussx[ix][iy][iz] = 0.0;
                gaussy[ix][iy][iz] = 0.0;
                gaussz[ix][iy][iz] = 0.0;
                gaussz_p[ix][iy][iz] = 0.0;
            }
        }
    }

    for(ipart=0; ipart<npart; ipart++) {
        ggem_sldq_to_glmesh(xb[ipart],fb[ipart],gaussx,gaussy,gaussz,gaussz_p,
                           weights[ipart],alpha,xz);
    }


    /* solve for the velocity and pressure at the mesh points */
    global_velocity_inhomogeneous(gaussx,gaussy,gaussz,gaussz_p,ubc1,ubc2,
                                  uxf,uyf,uzf,pf,
                                  dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
                                  xxi,cz, u1H,v1H,w1H,dwdz1H,dudz1H,dvdz1H,p1H,
                                  u2H,v2H,w2H,dwdz2H,dudz2H,dvdz2H,p2H);


    /* Interpolate velocity from the mesh to the bead position */
    for(ipart=0; ipart<npart; ipart++) {
        ggem_velocity_from_glmesh(xb[ipart],ub[ipart],uxf,uyf,uzf,xz);
    }

}

/******************************************************************************/

/* computes the traction at nodes due to the imposed flow:finf */
void ggem_calc_finf(double finf[beads][3],double xb[beads][3],
        double nrm[beads][3])
{
    int ibead;
    double pi;
    double U0;
    double gdot=GDOT;

    /* Compute the stress tensor and then traction at each node */
    if (SHEAR==0) {
        /* set U0=centerline velocity */
        U0 = gdot*Lz/4.0;
        for(ibead=0; ibead<beads; ibead++) {
            finf[ibead][0] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][0] + 4*U0/Lz*
                             (1-2*xb[ibead][2]/Lz)*nrm[ibead][2];
            finf[ibead][1] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][1];
            finf[ibead][2] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][2] + 4*U0/Lz*
                             (1-2*xb[ibead][2]/Lz)*nrm[ibead][0];
        }
    }
    else {
        for(ibead=0; ibead<beads; ibead++) {
            finf[ibead][0] = gdot*nrm[ibead][2];
            finf[ibead][1] = 0.0;
            finf[ibead][2] = gdot*nrm[ibead][0];
        }
    }
}

/******************************************************************************/

