#include "ggem.h"

/******************************************************************************/

void ggem_sldq_to_glmesh(GgemSolver* this, double* xb, double f[beads][3],
            double gaussx[nx][ny][nz], double gaussy[nx][ny][nz],
            double gaussz[nx][ny][nz], double gaussz_p[nx][ny][nz],
            double weight[beads], double alpha, double xz[nz])
{
    int ix, iy, iz;
    double dx, dy, dz;
    double xm[3],  r,  gr,  gr_p;
    double coeff;
    double cutoff;
    int ix0, iy0, iz0, iz1, iz2;
    int dx0, dy0;
    double rcut;
    int ixm, iym, izm;
    double exp_gr;
    double xzp1, xzm1;

    /* Mesh size */
    dx = Lx/nx;
    dy = Ly/ny;

    /* coeff of force density */
    static double coeff = pow(this->alpha/sqrt(M_PI), 3);

    /* find dx0, dy0, dz0 */
    static double rcut = 4.0/this->alpha;
    dx0 = rcut/dx;
    dy0 = rcut/dy;

    /* Distribute to mesh */
    int num_vetices = mesh_get_num_vertices(this->msh);

    for(int i=0; i < num_vertices; ++i) {
        ix0 = (int) (xb[3*i+0]/dx+0.5);
        iy0 = (int) (xb[3*i+1]/dy+0.5);

        /* top cutoff length */
        xzp1 = xb[3*i+2]+rcut;
        if (xzp1 >= Lz) {
            iz2 = nz-1;
        }
        else {
            iz2 = (int) (acos(1.0-2.0*xzp1/Lz)*(nz-1)/M_PI + 1.0); // +1 to be safe
        }

        /* bottom cutoff length */
        xzm1 = xb[3*i+2]-rcut;
        if (xzm1 <= 0.0) {
            iz1 = 0;
        }
        else {
            iz1 = (int) (acos(1.0-2.0*xzm1/Lz)*(nz-1)/M_PI - 1); // -1 to be safe
        }

        for(ix=ix0-dx0; ix<=ix0+dx0; ix++) {
            xm[0] = dx*ix - xb[3*i+0];

            for(iy=iy0-dy0; iy<=iy0+dy0; iy++) {
                xm[1] = dy*iy - xb[3*i+1];

                for(iz=iz1; iz<=iz2; iz++) {
                    xm[2] = xz[iz] - xb[3*i+2];

                    r = sqrt(xm[0]*xm[0]+xm[1]*xm[1]+xm[2]*xm[2]);

                    if ( r > rcut) {
                        continue;
                    }
                    exp_gr = exp(-alpha*alpha*r*r);
                    gr = coeff*exp_gr*(2.5-alpha*alpha*r*r);
                    gr_p = coeff*exp_gr*(3.5-alpha*alpha*r*r)*(-2.0*alpha*alpha*xm[2]);

                    /* Correct for periodicity */
                    if (ix >= nx) {
                        ixm = ix-nx;
                    }
                    else if (ix < 0) {
                        ixm = ix+nx;
                    }
                    else {
                        ixm = ix;
                    }

                    if (iy >= ny) {
                        iym = iy-ny;
                    }
                    else if (iy < 0) {
                        iym = iy+ny;
                    }
                    else {
                        iym = iy;
                    }

                    if (iz >= nz) {
                        continue;
                    }
                    else if (iz < 0) {
                        continue;
                    }
                    else {
                        izm = iz;
                    }

                    gaussx[ixm][iym][izm]   += gr*f[i][0]*weight[i];
                    gaussy[ixm][iym][izm]   += gr*f[i][1]*weight[i];
                    gaussz[ixm][iym][izm]   += gr*f[i][2]*weight[i];
                    gaussz_p[ixm][iym][izm] += gr_p*f[i][2]*weight[i];

                }
            }
        }
    }
}


/******************************************************************************/

/* Interpolate the global velocity at the field point  */
/* Uses quartic interpolation between the nodes of the containing cell */
/* currently assumes that the points are well inside the walls */
void ggem_velocity_from_glmesh(GgemSolver this, double* xb, double ub[beads][3],
            double uxf[nx][ny][nz], double uyf[nx][ny][nz],
            double uzf[nx][ny][nz], double xz[nz])
{
    int ix, iy, iz;
    int ix0, iy0, iz0, iz_check;
    int i, j, k;
    double dx, dy, dz;
    int ipoint;
    double ux1[5][5][5], uy1[5][5][5], uz1[5][5][5];
    double ux2[5][5], uy2[5][5], uz2[5][5];
    double ux3[5], uy3[5], uz3[5];
    double li[5];
    double ux, uy, uz;
    double x1, x2, x3, x4, x5, xc;
    int ibead;
    int kmap[5];
    double d_min=Lz,  dist;

    dx = Lx/nx;
    dy = Ly/ny;

    for(ibead=0; ibead<beads; ibead++) {

        ix0  = (int) (xb[ibead][0]/dx+0.5);
        iy0  = (int) (xb[ibead][1]/dy+0.5);
        iz0 = (int) (acos(1.0-2.0*xb[ibead][2]/Lz)*(nz-1)/pi);

        if (fabs(xz[iz0+1] - xb[ibead][2]) < fabs(xz[iz0] - xb[ibead][2])) {
            iz0 = iz0+1;
        }

        if (iz0 == 0) {
            kmap[0] = 0;
            kmap[1] = 1;
            kmap[2] = 2;
            kmap[3] = 3;
            kmap[4] = 4;
        }
        else if (iz0 == 1) {
            kmap[0] = -1;
            kmap[1] =  0;
            kmap[2] =  1;
            kmap[3] =  2;
            kmap[4] =  3;
        }
        else if (iz0 == nz-1) {
            kmap[0] = -4;
            kmap[1] = -3;
            kmap[2] = -2;
            kmap[3] = -1;
            kmap[4] =  0;
        }

        else if (iz0 == nz-2) {
            kmap[0] = -3;
            kmap[1] = -2;
            kmap[2] = -1;
            kmap[3] =  0;
            kmap[4] =  1;
        }
        else {
            kmap[0] = -2;
            kmap[1] = -1;
            kmap[2] = 0;
            kmap[3] = 1;
            kmap[4] = 2;
        }

        /* Store values of the velocity of the containing cell */
        for(i=-2; i<3; i++) {
            for(j=-2; j<3; j++) {
                for(k=-2; k<3; k++) {
                    ix = ix0+i;
                    iy = iy0+j;
                    iz = iz0+kmap[k+2];

                    /* correct for periodicity */
                    if (ix >= nx) {
                        ix = ix - nx;
                    }
                    else if (ix < 0) {
                        ix = ix + nx;
                    }

                    if (iy >= ny) {
                        iy = iy - ny;
                    }
                    else if (iy < 0) {
                        iy = iy + ny;
                    }

                    ux1[i+2][j+2][k+2]=uxf[ix][iy][iz];
                    uy1[i+2][j+2][k+2]=uyf[ix][iy][iz];
                    uz1[i+2][j+2][k+2]=uzf[ix][iy][iz];
                }
            }
        }

        /* Interpolate along x */
        x1 = (ix0-2)*dx;
        x2 = (ix0-1)*dx;
        x3 =  ix0*dx;
        x4 = (ix0+1)*dx;
        x5 = (ix0+2)*dx;
        xc = xb[ibead][0];

        // could remove reduntant calculation here--will do later
        li[0] = (xc-x2)*(xc-x3)*(xc-x4)*(xc-x5)/((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5));
        li[1] = (xc-x1)*(xc-x3)*(xc-x4)*(xc-x5)/((x2-x1)*(x2-x3)*(x2-x4)*(x2-x5));
        li[2] = (xc-x1)*(xc-x2)*(xc-x4)*(xc-x5)/((x3-x1)*(x3-x2)*(x3-x4)*(x3-x5));
        li[3] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x5)/((x4-x1)*(x4-x2)*(x4-x3)*(x4-x5));
        li[4] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x4)/((x5-x1)*(x5-x2)*(x5-x3)*(x5-x4));


        for(j=0; j<5; j++) {
            for(k=0; k<5; k++) {
                ux2[j][k] = 0.0;
                uy2[j][k] = 0.0;
                uz2[j][k] = 0.0;

                for(i=0; i<5; i++) {
                    ux2[j][k] += ux1[i][j][k]*li[i];
                    uy2[j][k] += uy1[i][j][k]*li[i];
                    uz2[j][k] += uz1[i][j][k]*li[i];
                }
            }
        }

        /* Interpolate along y */
        x1 = (iy0-2)*dy;
        x2 = (iy0-1)*dy;
        x3 =  iy0*dy;
        x4 = (iy0+1)*dy;
        x5 = (iy0+2)*dy;
        xc = xb[ibead][1];

        // could remove reduntant calculation here--will do later
        li[0] = (xc-x2)*(xc-x3)*(xc-x4)*(xc-x5)/((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5));
        li[1] = (xc-x1)*(xc-x3)*(xc-x4)*(xc-x5)/((x2-x1)*(x2-x3)*(x2-x4)*(x2-x5));
        li[2] = (xc-x1)*(xc-x2)*(xc-x4)*(xc-x5)/((x3-x1)*(x3-x2)*(x3-x4)*(x3-x5));
        li[3] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x5)/((x4-x1)*(x4-x2)*(x4-x3)*(x4-x5));
        li[4] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x4)/((x5-x1)*(x5-x2)*(x5-x3)*(x5-x4));


        for(k=0; k<5; k++) {
            ux3[k]=0.0;
            uy3[k]=0.0;
            uz3[k]=0.0;
            for(j=0; j<5; j++) {
                ux3[k] += ux2[j][k]*li[j];
                uy3[k] += uy2[j][k]*li[j];
                uz3[k] += uz2[j][k]*li[j];
            }
        }

        /* Interpolate along z */
        x1 = xz[(iz0+kmap[0])];
        x2 = xz[(iz0+kmap[1])];
        x3 = xz[(iz0+kmap[2])];
        x4 = xz[(iz0+kmap[3])];
        x5 = xz[(iz0+kmap[4])];
        xc = xb[ibead][2];

        li[0] = (xc-x2)*(xc-x3)*(xc-x4)*(xc-x5)/((x1-x2)*(x1-x3)*(x1-x4)*(x1-x5));
        li[1] = (xc-x1)*(xc-x3)*(xc-x4)*(xc-x5)/((x2-x1)*(x2-x3)*(x2-x4)*(x2-x5));
        li[2] = (xc-x1)*(xc-x2)*(xc-x4)*(xc-x5)/((x3-x1)*(x3-x2)*(x3-x4)*(x3-x5));
        li[3] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x5)/((x4-x1)*(x4-x2)*(x4-x3)*(x4-x5));
        li[4] = (xc-x1)*(xc-x2)*(xc-x3)*(xc-x4)/((x5-x1)*(x5-x2)*(x5-x3)*(x5-x4));

        ux = 0.0;
        uy = 0.0;
        uz = 0.0;

        for(k=0; k<5; k++) {
            ux += ux3[k]*li[k];
            uy += uy3[k]*li[k];
            uz += uz3[k]*li[k];
        }

        ub[ibead][0] += 8.0*M_PI*ux;
        ub[ibead][1] += 8.0*M_PI*uy;
        ub[ibead][2] += 8.0*M_PI*uz;
    }
}


