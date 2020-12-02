#include "ggem.h"
#include "utils_math.h"
#include <math.h>


/* positions and weights */
/* Non-singular */
double wphi_ns[NPHI],zphi_ns[NPHI];
double wr_ns[NR],zr_ns[NR];
/* Singular */
double wphi_s[NPHI],zphi_s[NPHI];
double wr_s[NR],zr_s[NR];

/******************************************************************************/

/* Calculates the area element */
double get_area_element(double x[3], double y[3])
{
    double z[3];

    cross(x, y, z);
    return norm(3, z);
}

/******************************************************************************/

/* Calculates the local Green function */
void calc_Gl(const double r[3], const double alpha, double* Gl)
{
    const int ncols = 3;
    double rg  = norm(3, r);
    double rg_sq = rg*rg;
    double rg_inv = 1.0/rg;
    double rg_inv_sq = rg_inv*rg_inv;
    double rg_inv_cu = rg_inv*rg_inv_sq;
    double erfc_rg = erfc(alpha*rg);
    double exp_rg = exp(-alpha*alpha*rg_sq);
    
    /* find the local green's function */
    Gl[ncols*0+0] = erfc_rg*(rg_inv+rg_inv_cu*r[0]*r[0])
                - M_2_SQRTPI*alpha*exp_rg*(1.0-rg_inv_sq*r[0]*r[0]);
    Gl[ncols*1+1] = erfc_rg*(rg_inv+rg_inv_cu*r[1]*r[1])
                - M_2_SQRTPI*alpha*exp_rg*(1.0-rg_inv_sq*r[1]*r[1]);
    Gl[ncols*2+2] = erfc_rg*(rg_inv+rg_inv_cu*xd[2]*r[2])
                - M_2_SQRTPI*alpha*exp_rg*(1.0-rg_inv_sq*r[2]*r[2]);
    
    Gl[ncols*0+1] = erfc_rg*(rg_inv_cu*r[0]*r[1])
                - M_2_SQRTPI*alpha*exp_rg*(-rg_inv_sq*r[0]*r[1]);
    Gl[ncols*0+2] = erfc_rg*(rg_inv_cu*r[0]*r[2])
                - M_2_SQRTPI*alpha*exp_rg*(-rg_inv_sq*r[0]*r[2]);
    Gl[ncols*1+2] = erfc_rg*(rg_inv_cu*r[1]*r[2])
                - M_2_SQRTPI*alpha*exp_rg*(-rg_inv_sq*r[1]*r[2]);
    
    Gl[ncols*1+0] = Gl[3*0+1];
    Gl[ncols*2+0] = Gl[3*0+2]
    Gl[ncols*2+1] = Gl[3*1+2];
}

/******************************************************************************/

void ggem_calc_sli_local(const double x0[3], const double x1[3],
        const double x2[3], const double x3[3], const double f1[3], 
        const double f2[3], const double f3[3], const double alpha, double u[3])
{
    int nphi=NPHI;
    int nr = NR;
    double r1,r2;
    double w;
    double phi,r;
    double l1,l2,l3;
    double xi;
    double eta;
    double omega;
    double x_xi[3];
    double x_eta[3];
    double f[3];
    double x[3];
    double r[3];
    double Gl[9];

    gauss_leg(nphi,wphi_ns,zphi_ns);
    gauss_leg(nr,wr_ns,zr_ns);

    /* find the area element omega */
    for(int i=0; i < 3; ++i) {
        x_xi[i] = x3[i] - x1[i];
        x_eta[i] = x2[i] - x1[i];
    }

    /* area element is constant for linear element */
    omega = get_area_element(x_xi, x_eta);

    /* initialize to zero */
    for(int i=0; i < 3; ++i) {
        u[i] = 0.0;
    }

    for(int iphi=0; iphi < nphi; ++iphi) {
        eta = 0.5*(1.0+zphi_ns[iphi]);

        for(int ir=0; ir < nr; ++ir) {
            xi = 0.5*(1.0+zr_ns[ir])*(1.0-eta);

            /* find the xi and eta coordinates */
            l1 = 1.0-xi-eta;
            l2 = eta;
            l3 = xi;

            /* coefficient of green's function at the current point */
            for(int i=0; i<3; ++i) {
                f[i] = l1*f1[i] + l2*f2[i] + l3*f3[i];
            }

            /* position of current point */
            for(int i=0; i < 3; ++i) {
                x[i] = l1*x1[i] + l2*x2[i] + l3*x3[i];
            }

            /* separation vector: field - pole */
            for(int i=0; i < 3; ++i) {
                r[i] = x0[i] - x[i];
            }

            /* Correct for periodicity */
            r[0] = r[0] - Lx*floor(r[0]/Lx+0.5);
            r[1] = r[1] - Ly*floor(r[1]/Ly+0.5);

            /* global coordinates distance rg */
            calc_Gl(r, alpha, Gl);
            /* weight */
            weight  = wphi_ns[iphi]*wr_ns[ir]*(1/2.0)*(1-eta)/2.0*omega;

            /* find contribution to velocity */
            for(int i=0; i < 3; ++i) {
                for(int j=0; j < 3; ++j) {
                    u[i] += f[j]*Gl[i][j]*weight;
                }
            }

        }
    }
}

/******************************************************************************/

void ggem_calc_sli_local_polar(const double x0[3], const double x1[3],
        const double x2[3], const double x3[3], const double f1[3], 
        const double f2[3], const double f3[3], const double alpha, double u[3])
{
    int nphi=NPHI;
    int nr = NR;
    double phi1, phi2;
    double r1, r2;
    double w;
    double phi, r;
    double l1, l2, l3;
    double xi;
    double eta;
    double f[3];
    double x[3];
    double omega;
    double x_xi[3];
    double x_eta[3];
    double r[3];
    double Gl[9];

    phi1 = 0.0;
    phi2 = M_PI_2;

    /* find the area element omega */
    for(int i=0; i<3; ++i) {
        x_xi[i]  = x3[i] - x1[i];
        x_eta[i] = x2[i] - x1[i];
    }
    // area element is constant for linear element
    omega = get_area_element(x_xi, x_eta); 

    /* initialize to zero */
    for(int i=0; i<3; ++i) {
        u[i]=0.0;
    }

    for(int iphi=0; iphi<nphi; ++iphi) {
        phi = phi1 + 0.5*(1.0+zphi_s[iphi])*(phi2-phi1);
        r1 = 0.0;
        r2 = sin(M_PI/4.0)/sin(3*M_PI/4.0-phi);

        for(int ir = 0; ir < nr; ++ir) {
            r = r1 + 0.5*(1.0+zr_s[ir])*(r2-r1);

            /* find the xi and eta coordinates */
            xi = r*sin(phi);
            eta = r*cos(phi);

            /* find the xi and eta coordinates */
            l1 = 1.0-xi-eta;
            l2 = eta;
            l3 = xi;

            /* coefficient of green's function at the current point */
            for(int i=0; i<3; ++i) {
                f[i] = l1*f1[i] + l2*f2[i] + l3*f3[i];
            }

            /* position of current point */
            for(int i=0; i<3; ++i) {
                x[i] = l1*x1[i]+ l2*x2[i] + l3*x3[i];
            }

            /* separation vector: field - pole */
            for(int i=0; i<3; ++i) {
                r[i] = x0[i]-x[i];
            }

            /* Correct for periodicity */
            r[0] = r[0] - Lx*floor(xd[0]/Lx+0.5);
            r[1] = r[1] - Ly*floor(xd[1]/Ly+0.5);

            calc_Gl(r, alpha, Gl);

            /* weight */
            w  = wphi_s[iphi]*wr_s[ir]*(pi/2.0/2.0)*(r2-r1)/2.0*omega;

            /* find contribution to velocity */
            for(int i=0; i<3; ++i) {
                for(int j=0; j<3; ++j) {
                    u[i] += f[j]*Gl[3*i+j]*w*r;
                }
            }
        }
    }
}

/******************************************************************************/
