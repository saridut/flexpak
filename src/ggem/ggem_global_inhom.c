/* Solves the velocity field for the global force density */

#include "ggem.h"
#include <stdlib.h>
#include <stdio.h>

#define LAMBDA 0
#define STRESS_HI 0

/* solves the system: a11y" + a22y = f */
/* algorithm given in appexdix B of Peyret's book */
/* Also scales appropriately for mapping from (1,-1) to (0,Lz) */

void solve_quasi_tridiagonal(double nu, double b,double complex fhat[nz],
            double complex cheb_bc[2], double complex uhat[nz], double cz[nz])
{
    double e[nz+4];
    double p[nz],q[nz],r[nz];
    int Ie,Io;
    int i,k;
    double complex f[nz];
    double complex X[nz],Y[nz];
    double complex theta[nz+1],lambda[nz+1];
    double complex sum1,sum2;
    double complex g;
    double complex w[nz];

    nu = nu*4.0/Lz/Lz; // account for mapping double derivative.

    /* set e (in global coordinates) */
    for(i=0; i<nz; i++) {
        e[i]=1.0;
    }
    e[nz]=0.0;
    e[nz+1]=0.0;
    e[nz+2]=0.0;
    e[nz+3]=0.0;

    /* set p[k],r[k],q[k]  (in global coordinates)*/
    for(i=2; i<nz; i++) {
        p[i] = cz[i-2]/4.0/i/(i-1);
        q[i] = -e[i+2]/2.0/(i*i-1);
        r[i] = e[i+4]/4.0/i/(i+1);
    }

    Ie = nz - nz/2;  // # of even equations
    Io = nz - Ie;    // # of odd equations

    Ie = Ie - 1; // subtracting eqs for the BC
    Io = Io - 1; // subtracting eqs for the BC


    /* --------- solve the even set of equations equation ----- */

    /* set the rhs of the tridiagonal solve */
    for(i=0; i<Ie-1; i++) { // first Ie-1 equations
        k = 2+2*i; // global eq #
        f[i] = p[k]*fhat[k-2] + q[k]*fhat[k] + r[k]*fhat[k+2];
    }
    i = Ie-1;
    k = 2+2*i;
    f[i] =  p[k]*fhat[k-2] + q[k]*fhat[k];

    /* find X & Y */
    i = Ie-1;
    k = 2+2*i;
    X[Ie-1] = -p[k]*b/(b*q[k]-nu);
    Y[Ie-1] = f[i]/(b*q[k]-nu);

    for(i=Ie-2; i>=0; i--) {
        k = 2+2*i;
        X[i] = -p[k]*b/(q[k]*b-nu + b*r[k]*X[i+1]);
        Y[i] = (f[i]-b*r[k]*Y[i+1])/(q[k]*b-nu + b*r[k]*X[i+1]);
    }

    /* find theta, lambda */
    theta[0]=1;
    lambda[0]=0;
    for(i=1; i<=Ie; i++) {
        theta[i] = X[i-1]*theta[i-1];
        lambda[i] = X[i-1]*lambda[i-1] + Y[i-1];
    }

    /* calculate w0 */
    sum1=0.0;
    for(i=0; i<=Ie; i++) {
        sum1 += theta[i];
    }
    sum2=0.0;
    for(i=1; i<=Ie; i++) {
        sum2 += lambda[i];
    }

    g = (cheb_bc[0]+cheb_bc[1])/2.0;
    w[0] = (g - sum2)/sum1;

    for(i=1; i<=Ie; i++) {
        w[i] = X[i-1]*w[i-1] + Y[i-1];
    }

    /* store in uhat */
    for(i=0; i<=Ie; i++) {
        k = 2*i;
        uhat[k] = w[i];
    }


    /* --------- solve the odd set of equations equation ----- */

    /* set the rhs of the tridiagonal solve */
    for(i=0; i<Io-1; i++) { // first Ie-1 equations
        k = 3+2*i; // global eq #
        f[i] = p[k]*fhat[k-2] + q[k]*fhat[k] + r[k]*fhat[k+2];
    }
    i = Io-1;
    k = 3+2*i;
    f[i] =  p[k]*fhat[k-2] + q[k]*fhat[k];

    /* find X & Y */
    i = Io-1;
    k = 3+2*i;
    X[Io-1] = -p[k]*b/(b*q[k]-nu);
    Y[Io-1] = f[i]/(b*q[k]-nu);

    for(i=Io-2; i>=0; i--) {
        k = 3+2*i;
        X[i] = -p[k]*b/(q[k]*b-nu + b*r[k]*X[i+1]);
        Y[i] = (f[i]-b*r[k]*Y[i+1])/(q[k]*b-nu + b*r[k]*X[i+1]);
    }

    /* find theta, lambda */
    theta[0]=1;
    lambda[0]=0;
    for(i=1; i<=Io; i++) {
        theta[i] = X[i-1]*theta[i-1];
        lambda[i] = X[i-1]*lambda[i-1] + Y[i-1];
    }

    /* calculate w0 */
    sum1=0.0;
    for(i=0; i<=Io; i++) {
        sum1 += theta[i];
    }
    sum2=0.0;
    for(i=1; i<=Io; i++) {
        sum2 += lambda[i];
    }

    g = (cheb_bc[1]-cheb_bc[0])/2.0;

    w[0] = (g - sum2)/sum1;
    for(i=1; i<=Io; i++) {
        w[i] = X[i-1]*w[i-1] + Y[i-1];
    }

    /* store in uhat */
    for(i=0; i<=Io; i++) {
        k = 2*i+1;
        uhat[k] = w[i];
    }

}

/*-------------------------------------------------------------
 * Computes the first derivative p_i(1) in terms of p_k
 * Also scales appropriately for mapping from (1,-1) to 0,Lz
 * -----------------------------------------------------------*/
void compute_coef_derivative(double complex u[nz],double complex uprime[nz],
                             double cz[nz])
{
    int i;
    uprime[nz-1]=0.0;

    uprime[nz-2]=2.0*(nz-1)*u[nz-1]/cz[nz-2];

    for(i=nz-3; i>=0; i--) {
        uprime[i] = (uprime[i+2] + 2.0*(i+1)*u[i+1])/cz[i];
    }


    /* scale */
    for(i=0; i<nz; i++) {
        uprime[i] = uprime[i]*(-2.0/Lz);
    }

}

/*------------------------------------------------------------------------------------
 -------------- Solve for the global velocity and pressure -------------------------
 -----------------------------------------------------------------------------------*/

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
        double complex ptr2[nx][ny/2+1][nz])
{


    fftw_plan plan;
    double k1,k2,ksq;
    double tmp1[nx][ny];
    double complex tmp2[nx][ny/2+1];
    double complex bcx1[nx][ny/2+1], bcy1[nx][ny/2+1], bcz1[nx][ny/2+1];
    double complex bcx2[nx][ny/2+1], bcy2[nx][ny/2+1], bcz2[nx][ny/2+1];
    double complex rhox[nx][ny/2+1][nz], rhoy[nx][ny/2+1][nz],
           rhoz[nx][ny/2+1][nz], rhoz_p[nx][ny/2+1][nz];
    double complex ptr[nx][ny/2+1][nz];
    double complex phat[nz],pdhat[nz],uhat[nz],vhat[nz],what[nz],wdhat[nz],
           udhat[nz],vdhat[nz];
    double complex rhox_hat[nz],rhoy_hat[nz],rhoz_hat[nz],rhoz_p_hat[nz];
    double complex pbc1[nx][ny/2+1],pbc2[nx][ny/2+1]; // pressure BC
    double complex A[2][2],b[2];
    double complex pb1,pb2;
    double complex dudz_t[nx][ny/2+1][nz],dvdz_t[nx][ny/2+1][nz];


    int i,j,k,n1,n2;
    int info;
    double pi;
    double dz,dzsq,KX,KY,eta=1;
    FILE *fp;
    int il;

    double complex fhat[nz];
    double complex cheb_in[2*nz-2], cheb_out[2*nz-2];
    int Nz=2*nz-2; // Number of Fourier points for computing chebyshev FFT
    double complex cheb_bc[2];
    double a11,a22;

    double complex p2dhat[nz],p2d[nz],rhs[nz];

    pi = 4.0*atan(1.0);

    dz = Lz/(nz-1);
    dzsq = dz*dz;
    KX = 2*pi/Lx;
    KY = 2*pi/Ly;


    /*-------------------------------------------------------------*/
    /* Transform the u_far boundary conditions */
    /*-------------------------------------------------------------*/

    /* Bottom wall boundary conditions */
    /* x component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc1[i][j][0];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcx1,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* y component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc1[i][j][1];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcy1,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* z component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc1[i][j][2];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcz1,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* pressure BC == 0, no need to transform */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            pbc1[n1][n2]=0.0;
        }
    }



    /* Top wall boundary conditions */
    /* x component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc2[i][j][0];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcx2,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* y component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc2[i][j][1];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcy2,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* z component of velocity */
    for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
            tmp1[i][j] = ubc2[i][j][2];
        }
    }
    plan = fftw_plan_dft_r2c_2d(nx,ny,tmp1,bcz2,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* pressure BC == 0, no need to transform */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            pbc2[n1][n2]=0.0;
        }
    }

    /*-------------------------------------------------------------*/
    /* Transform the Gaussian force distribution to fourier-space */
    /*-------------------------------------------------------------*/

    for(k=0; k<nz; k++) {
        /* transform x component of density */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                tmp1[i][j] = gaussx[i][j][k];
            }
        }

        plan = fftw_plan_dft_r2c_2d(nx, ny, tmp1, tmp2, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store in rhox */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhox[n1][n2][k] = tmp2[n1][n2];
            }
        }

        /* transform y component of density */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                tmp1[i][j] = gaussy[i][j][k];
            }
        }

        plan= fftw_plan_dft_r2c_2d(nx, ny, tmp1, tmp2, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store in rhoy */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoy[n1][n2][k] = tmp2[n1][n2];
            }
        }

        /* transform z component of density */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                tmp1[i][j] = gaussz[i][j][k];
            }
        }

        plan= fftw_plan_dft_r2c_2d(nx, ny, tmp1, tmp2, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store in rhoz */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoz[n1][n2][k] = tmp2[n1][n2];
            }
        }

        /* transform the z (y in text) derivative of the rhoz */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                tmp1[i][j] = gaussz_p[i][j][k];
            }
        }

        plan= fftw_plan_dft_r2c_2d(nx, ny, tmp1, tmp2, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store in rhoz */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoz_p[n1][n2][k] = tmp2[n1][n2];
            }
        }
    }




    /* --------------------------------------------------------------------*/
    /*------------ Solution of the equation in Z for each fourier mode */
    /*------------------------------------------------------------------*/
    for(n1=0; n1<nx; n1++) {
        if (n1<=nx/2) {
            k1 = n1;
        }
        else {
            k1 = n1 - nx;
        }

        k1 = k1*KX;

        /* Real data, so n/2+1 terms needed, remaining frequencies are automatically obtained
        from hermitian symmetry */
        for(n2=0; n2<ny/2+1; n2++) {
            k2 = n2*KY;

            // don't forget to scale k1 and k2
            ksq = k1*k1 + k2*k2;   //k^2;



            /*---------------------------------------------------------------
             ------------- Preliminaries for all HELMHOLTZ EQUATIONs---------
             ----------------------------------------------------------------*/

            /* -----find chebyshev coefficient of rhox-------- */
            for(i=0; i<nz; i++) {
                cheb_in[i] = rhox[n1][n2][i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhox_hat[i] = cheb_out[i]/(nz-1)/cz[i];
            }


            /* -----find chebyshev coefficient of rhoy-------- */
            for(i=0; i<nz; i++) {
                cheb_in[i] = rhoy[n1][n2][i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoy_hat[i] = cheb_out[i]/(nz-1)/cz[i];
            }

            /* -----find chebyshev coefficient of rhoz-------- */
            for(i=0; i<nz; i++) {
                cheb_in[i] = rhoz[n1][n2][i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoz_hat[i] = cheb_out[i]/(nz-1)/cz[i];
            }

            /* -----find chebyshev coefficient of rhoz_p-------- */
            for(i=0; i<nz; i++) {
                cheb_in[i] = rhoz_p[n1][n2][i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoz_p_hat[i] = cheb_out[i]/(nz-1)/cz[i];
            }


            /*---------------------------------------------------------------
             ------------------- PRESSURE HELMHOLTZ EQUATION ----------------
             ----------------------------------------------------------------*/
            /* Find the rhs of the helmholtz equation */
            for(i=0; i<nz; i++) {
                fhat[i] = k1*rhox_hat[i]*I + k2*rhoy_hat[i]*I + rhoz_p_hat[i];
            }

            /* set the bc before solving */
            cheb_bc[0]=pbc1[n1][n2];
            cheb_bc[1]=pbc2[n1][n2];

            /* set coefficient of the helmholtz equation */
            a11 = -1.0;
            a22 = -ksq;
            solve_quasi_tridiagonal(a11,a22,fhat,cheb_bc,phat,cz);

            /* compute chebyshev coefficient of derivative of pressure */
            compute_coef_derivative(phat,pdhat,cz);

            /*---- compute fourier coefficients of pressure at mesh points from phat using FFT ---*/
            for(i=0; i<nz; i++) {
                cheb_in[i] = phat[i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* scale 0 and Nth by 2.0*/
            cheb_in[0]    = 2.0*cheb_in[0];
            cheb_in[nz-1] = 2.0*cheb_in[nz-1];

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                ptr[n1][n2][i] = cheb_out[i]/2.0;
            }


            /*---------------------------------------------------------------
             ------------------- U HELMHOLTZ EQUATION ----------------------
             ----------------------------------------------------------------*/
            /* Find the rhs of the helmholtz equation */
            for(i=0; i<nz; i++) {
                fhat[i] = -k1*phat[i]*I + rhox_hat[i];
            }

            /* set the bc before solving */
            cheb_bc[0]=bcx1[n1][n2];
            cheb_bc[1]=bcx2[n1][n2];

            /* set coefficient of the helmholtz equation */
            a11 = 1.0;
            a22 = ksq;
            solve_quasi_tridiagonal(a11,a22,fhat,cheb_bc,uhat,cz);

#if (LAMBDA==0 && STRESS_HI == 1)
            /* compute chebyshev coefficient of derivative of u */
            compute_coef_derivative(uhat,udhat,cz);
#endif

            /*---- compute fourier coefficients of velocity at mesh points from uhat using FFT ---*/
            for(i=0; i<nz; i++) {
                cheb_in[i] = uhat[i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* scale 0 and Nth by 2.0*/
            cheb_in[0]    = 2.0*cheb_in[0];
            cheb_in[nz-1] = 2.0*cheb_in[nz-1];

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhox[n1][n2][i] = cheb_out[i]/2.0;
            }

#if (LAMBDA==0 && STRESS_HI == 1)
            /*---- compute fourier coefficients of velocity derivative at mesh points from udhat using FFT ---*/
            for(i=0; i<nz; i++) {
                cheb_in[i] = udhat[i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* scale 0 and Nth by 2.0*/
            cheb_in[0]    = 2.0*cheb_in[0];
            cheb_in[nz-1] = 2.0*cheb_in[nz-1];

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                dudz_t[n1][n2][i] = cheb_out[i]/2.0;
            }
#endif

            /*---------------------------------------------------------------
             ------------------- V HELMHOLTZ EQUATION ----------------------
             ----------------------------------------------------------------*/
            /* Find the rhs of the helmholtz equation */
            for(i=0; i<nz; i++) {
                fhat[i] = -k2*phat[i]*I + rhoy_hat[i];
            }

            /* set the bc before solving */
            cheb_bc[0]=bcy1[n1][n2];
            cheb_bc[1]=bcy2[n1][n2];

            /* set coefficient of the helmholtz equation */
            a11 = 1.0;
            a22 = ksq;
            solve_quasi_tridiagonal(a11,a22,fhat,cheb_bc,vhat,cz);

#if (LAMBDA==0 && STRESS_HI == 1)
            /* compute chebyshev coefficient of derivative of v */
            compute_coef_derivative(vhat,vdhat,cz);
#endif

            /*---- compute fourier coefficients of velocity at mesh points from vhat using FFT ---*/
            for(i=0; i<nz; i++) {
                cheb_in[i] = vhat[i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* scale o and Nth by 2.0*/
            cheb_in[0]    = 2.0*cheb_in[0];
            cheb_in[nz-1] = 2.0*cheb_in[nz-1];

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoy[n1][n2][i] = cheb_out[i]/2.0;
            }

#if (LAMBDA==0 && STRESS_HI == 1)
            /*---- compute fourier coefficients of velocity derivative at mesh points from vdhat using FFT ---*/
            for(i=0; i<nz; i++) {
                cheb_in[i] = vdhat[i];
            }
            for(i=nz; i<2*nz-2; i++) {
                cheb_in[i] = cheb_in[2*nz-2-i];
            }

            /* scale o and Nth by 2.0*/
            cheb_in[0]    = 2.0*cheb_in[0];
            cheb_in[nz-1] = 2.0*cheb_in[nz-1];

            /* Perform FFT */
            plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                    FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);

            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                dvdz_t[n1][n2][i] = cheb_out[i]/2.0;
            }
#endif
            /*---------------------------------------------------------------
            	 ------------------- W HELMHOLTZ EQUATION ----------------------
            	 ----------------------------------------------------------------*/
            if (n1 > 0 || n2 > 0) {
                /* Find the rhs of the helmholtz equation */
                for(i=0; i<nz; i++) {
                    fhat[i] = -pdhat[i] + rhoz_hat[i];
                }

                /* set the bc before solving */
                cheb_bc[0]=bcz1[n1][n2];
                cheb_bc[1]=bcz2[n1][n2];

                /* set coefficient of the helmholtz equation */
                a11 = 1.0;
                a22 = ksq;
                solve_quasi_tridiagonal(a11,a22,fhat,cheb_bc,what,cz);

                /* compute chebyshev coefficient of derivative of w */
                compute_coef_derivative(what,wdhat,cz);

                /*---- compute fourier coefficients of velocity at mesh points from what using FFT ---*/
                for(i=0; i<nz; i++) {
                    cheb_in[i] = what[i];
                }
                for(i=nz; i<2*nz-2; i++) {
                    cheb_in[i] = cheb_in[2*nz-2-i];
                }

                /* scale 0 and Nth by 2.0*/
                cheb_in[0]    = 2.0*cheb_in[0];
                cheb_in[nz-1] = 2.0*cheb_in[nz-1];

                /* Perform FFT */
                plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                        FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);

                /* scale and store the relevant first nz components */
                for(i=0; i<nz; i++) {
                    rhoz[n1][n2][i] = cheb_out[i]/2.0;
                }

                /*---- compute fourier coefficients of velocity derivative at mesh points from uhat using FFT ---*/
                for(i=0; i<nz; i++) {
                    cheb_in[i] = wdhat[i];
                }
                for(i=nz; i<2*nz-2; i++) {
                    cheb_in[i] = cheb_in[2*nz-2-i];
                }

                /* scale o and Nth by 2.0*/
                cheb_in[0]    = 2.0*cheb_in[0];
                cheb_in[nz-1] = 2.0*cheb_in[nz-1];

                /* Perform FFT */
                plan = fftw_plan_dft_1d(2*nz-2, cheb_in, cheb_out, FFTW_FORWARD,
                                        FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);

                /* scale and store the relevant first nz components */
                for(i=0; i<nz; i++) {
                    rhoz_p[n1][n2][i] = cheb_out[i]/2.0;
                }

            }
            else { // k == 0 term, vg cancels vl => set it to equal the BC */

                for(i=0; i<nz; i++) {
                    rhoz[n1][n2][i] = 0.5*(bcz1[n1][n2]+bcz2[n1][n2]);
                }

                for(i=0; i<nz; i++) {
                    rhoz_p[n1][n2][i] = 0.0;
                }
            }


            /* solve the influence matrix */
            if (n1 > 0 || n2 > 0) {
                A[0][0] = k1*rhox1[n1][n2][0]*I    + k2*rhoy1[n1][n2][0]*I    +
                          rhoz_p1[n1][n2][0];
                A[1][0] = k1*rhox1[n1][n2][nz-1]*I + k2*rhoy1[n1][n2][nz-1]*I +
                          rhoz_p1[n1][n2][nz-1];

                A[0][1] = k1*rhox2[n1][n2][0]*I    + k2*rhoy2[n1][n2][0]*I    +
                          rhoz_p2[n1][n2][0];
                A[1][1] = k1*rhox2[n1][n2][nz-1]*I + k2*rhoy2[n1][n2][nz-1]*I +
                          rhoz_p2[n1][n2][nz-1];

                b[0] = -(k1*rhox[n1][n2][0]*I    + k2*rhoy[n1][n2][0]*I    +
                         rhoz_p[n1][n2][0]);
                b[1] = -(k1*rhox[n1][n2][nz-1]*I + k2*rhoy[n1][n2][nz-1]*I + rhoz_p[n1][n2][nz
                         -1]);

                pb1 = (A[0][1]*b[1] - A[1][1]*b[0])/(A[0][1]*A[1][0]-A[0][0]*A[1][1]);
                pb2 = (A[1][0]*b[0] - A[0][0]*b[1])/(A[0][1]*A[1][0]-A[0][0]*A[1][1]);
            }
            else { // this just sets the isotropic pressure--set to zero
                pb1=1.0;
                pb2=1.0;
            }
            //pb1=0.0;pb2=0.0;
            //printf("%E\t%E\n",creal(b[0]),creal(b[1]));

            /* find the overall solution */
            //printf("Before:%d\t%d\t%E\t%E\n",n1,n2,creal(rhox[0][1][0]),cimag(rhox[0][1][0]);
            //printf("After:%d\t%d\t%E\n",n1,n2,creal(rhox[0][1][0]));

            //	printf("Before:%d\t%d\t%E\n",n1,n2,creal(rhox[n1][n2][n2/2]));
            for(i=0; i<nz; i++) {

                rhox[n1][n2][i] = rhox[n1][n2][i] + pb1*rhox1[n1][n2][i] +
                                  pb2*rhox2[n1][n2][i];
                rhoy[n1][n2][i] = rhoy[n1][n2][i] + pb1*rhoy1[n1][n2][i] +
                                  pb2*rhoy2[n1][n2][i];
                rhoz[n1][n2][i] = rhoz[n1][n2][i] + pb1*rhoz1[n1][n2][i] +
                                  pb2*rhoz2[n1][n2][i];
                ptr[n1][n2][i]  = ptr[n1][n2][i]  + pb1*ptr1[n1][n2][i]  + pb2*ptr2[n1][n2][i];
#if (LAMBDA==0 && STRESS_HI == 1)
                rhoz_p[n1][n2][i] = rhoz_p[n1][n2][i] + pb1*rhoz_p1[n1][n2][i] +
                                    pb2*rhoz_p2[n1][n2][i];
                dudz_t[n1][n2][i] = dudz_t[n1][n2][i] + pb1*dudz1[n1][n2][i] +
                                    pb2*dudz2[n1][n2][i];
                dvdz_t[n1][n2][i] = dvdz_t[n1][n2][i] + pb1*dvdz1[n1][n2][i] +
                                    pb2*dvdz2[n1][n2][i];
#endif
            }
            //	printf("After:%d\t%d\t%E\n",n1,n2,creal(rhox[n1][n2][n2/2]));

        }
    }

    /*-------------------------------------------------------------*/
    /* Solve for velocity on the mesh using inverse transform */
    /*-------------------------------------------------------------*/
    /* At each point along the height of the slit find the velocity */
    for(k=0; k<nz; k++) {

        /*---------- ux in real space ----------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = rhox[n1][n2][k];
            }
        }


        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                uxf[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* -------------- uy in real space --------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = rhoy[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                uyf[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }


        /* ------------- uz in real space---------------------*/
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = rhoz[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                uzf[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ------------------ pressure in real space----------------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = ptr[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                pf[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

#if (LAMBDA==0 && STRESS_HI == 1)
        /* ----------------------- dudx in real space ------------------------ */
        for(n1=0; n1<nx; n1++) {
            if (n1<=nx/2) {
                k1 = n1;
            }
            else {
                k1 = n1 - nx;
            }
            k1 = k1*KX;

            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = k1*rhox[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dudx[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ----------------- dudy in real space ----------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                k2 = n2*KY;
                tmp2[n1][n2] = k2*rhox[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dudy[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ---------------- dudz in real space ---------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = dudz_t[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dudz[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ----------------------- dvdx in real space ------------------------ */
        for(n1=0; n1<nx; n1++) {
            if (n1<=nx/2) {
                k1 = n1;
            }
            else {
                k1 = n1 - nx;
            }
            k1 = k1*KX;

            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = k1*rhoy[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dvdx[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ----------------- dvdy in real space ----------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                k2 = n2*KY;
                tmp2[n1][n2] = k2*rhoy[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dvdy[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ---------------- dvdz in real space ---------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = dvdz_t[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dvdz[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ----------------------- dwdx in real space ------------------------ */
        for(n1=0; n1<nx; n1++) {
            if (n1<=nx/2) {
                k1 = n1;
            }
            else {
                k1 = n1 - nx;
            }
            k1 = k1*KX;

            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = k1*rhoz[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dwdx[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ----------------- dwdy in real space ----------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                k2 = n2*KY;
                tmp2[n1][n2] = k2*rhoz[n1][n2][k]*I;
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dwdy[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }

        /* ---------------- dwdz in real space ---------------- */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                tmp2[n1][n2] = rhoz_p[n1][n2][k];
            }
        }

        /* peform transform */
        plan=fftw_plan_dft_c2r_2d(nx,ny,tmp2,tmp1,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        /* store results while scaling to account for the lack of normalization in fftw */
        for(i=0; i<nx; i++) {
            for(j=0; j<ny; j++) {
                dwdz[i][j][k] = tmp1[i][j]/nx/ny;
            }
        }
#endif


    }



}
