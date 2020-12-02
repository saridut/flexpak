#include "ggem.h"
#include <stdlib.h>
#include <stdio.h>

/******************************************************************************/
/* Solves the velocity field for the global force density */

void ggem_global_velocity_homogeneous(
        double complex rhox[nx][ny/2+1][nz], 
        double complex rhoy[nx][ny/2+1][nz],
        double complex rhoz[nx][ny/2+1][nz],
        double complex rhoz_p[nx][ny/2+1][nz],
	    double complex dudz[nx][ny/2+1][nz],
        double complex dvdz[nx][ny/2+1][nz],
	    double complex ptr[nx][ny/2+1][nz],
        double cz[nz],
        int problem_no)
{
    fftw_plan plan;
    double k1, k2, ksq;
    double tmp1[nx][ny];

    double complex tmp2[nx][ny/2+1];

    double complex bcx1[nx][ny/2+1];
    double complex bcy1[nx][ny/2+1];

    double complex bcz1[nx][ny/2+1];
    double complex bcx2[nx][ny/2+1];
    double complex bcy2[nx][ny/2+1];
    double complex bcz2[nx][ny/2+1];

    double complex phat[nz];
    double complex pdhat[nz];
    double complex uhat[nz];
    double complex vhat[nz];
    double complex what[nz];
    double complex wdhat[nz]; 
    double complex udhat[nz];
    double complex vdhat[nz];

    double complex rhox_hat[nz];
    double complex rhoy_hat[nz];
    double complex rhoz_hat[nz];
    double complex rhoz_p_hat[nz];

    double complex pbc1[nx][ny/2+1]; //pressure BC
    double complex pbc2[nx][ny/2+1]; // pressure BC

    int i, j, k, n1, n2;
    int info;
    double pi;
    double dz, dzsq, KX, KY, eta=1;
    int il;

    double complex fhat[nz];
    double complex cheb_in[2*nz-2], cheb_out[2*nz-2];
    int Nz=2*nz-2; // Number of Fourier points for computing chebyshev FFT
    double complex cheb_bc[2];
    double a11, a22;

    pi = M_PI;
    dz = Lz/(nz-1);
    dzsq = dz*dz;
    KX = 2*pi/Lx;
    KY = 2*pi/Ly;

    /* Bottom wall boundary conditions */
    /* x component of velocity == 0 */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcx1[n1][n2] = 0.0;
        }
    }

    /* y component of velocity */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcy1[n1][n2] = 0.0;
        }
    }

    /* z component of velocity */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcz1[n1][n2] = 0.0;
        }
    }

    /* pressure BC == 0 or 1 no need to transform */
    if (problem_no == 1) {
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                pbc1[n1][n2]=1.0;
            }
        }
    }
    else {
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                pbc1[n1][n2]=0.0;
            }
        }
    }

    /* Top wall boundary conditions */
    /* x component of velocity */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcx2[n1][n2] = 0.0;
        }
    }

    /* y component of velocity */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcy2[n1][n2] = 0.0;
        }
    }

    /* z component of velocity */
    for(n1=0; n1<nx; n1++) {
        for(n2=0; n2<ny/2+1; n2++) {
            bcz2[n1][n2] = 0.0;
        }
    }

    /* pressure BC == 0 or 1 no need to transform */
    if (problem_no == 2) {
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                pbc2[n1][n2]=1.0;
            }
        }
    }
    else {

        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                pbc2[n1][n2]=0.0;
            }
        }
    }

    /*-------------------------------------------------------------*/
    /* Transform the Gaussian force distribution to fourier-space */
    /*-------------------------------------------------------------*/
    for(k=0; k<nz; k++) {
        /* transform x component of density == 0*/
        /* store in rhox */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhox[n1][n2][k] = 0.0;
            }
        }

        /* transform y component of density */
        /* store in rhoy */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoy[n1][n2][k] = 0.0;
            }
        }

        /* transform z component of density */
        /* store in rhoz */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoz[n1][n2][k] = 0.0;
            }
        }

        /* transform the z (y in text) derivative of the rhoz */
        /* store in rhoz */
        for(n1=0; n1<nx; n1++) {
            for(n2=0; n2<ny/2+1; n2++) {
                rhoz_p[n1][n2][k] = 0.0;
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
            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhox_hat[i] = 0.0;
            }

            /* -----find chebyshev coefficient of rhoy-------- */
            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoy_hat[i] = 0.0;
            }

            /* -----find chebyshev coefficient of rhoz-------- */
            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoz_hat[i] = 0.0;
            }

            /* -----find chebyshev coefficient of rhoz_p-------- */
            /* scale and store the relevant first nz components */
            for(i=0; i<nz; i++) {
                rhoz_p_hat[i] = 0.0;
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

            /* find chebyshev coefficient of rhoz_p */
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
            solve_quasi_tridiagonal(a11, a22, fhat, cheb_bc, uhat, cz);

            /* compute chebyshev coefficient of derivative of u */
            compute_coef_derivative(uhat, udhat, cz);

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

            /*---- compute fourier coefficients of velocity derivative at mesh points from uhat using FFT ---*/
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
                dudz[n1][n2][i] = cheb_out[i]/2.0;
            }

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

            /* compute chebyshev coefficient of derivative of u */
            compute_coef_derivative(vhat,vdhat,cz);

            /*---- compute fourier coefficients of velocity at mesh points from uhat using FFT ---*/
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

            /*---- compute fourier coefficients of velocity derivative at mesh points from vhat using FFT ---*/
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
                dvdz[n1][n2][i] = cheb_out[i]/2.0;
            }

            /*---------------------------------------------------------------
              ------------------- W HELMHOLTZ EQUATION ----------------------
              ----------------------------------------------------------------*/
            if (1) { //(n1 > 0 || n2 > 0) //k=0 doesn't satisfy continuity; + and - solution must have equal coeff for k=0
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

                //compute fourier coefficients of velocity at mesh points from
                //uhat using FFT
                for(i=0; i<nz; i++) {
                    cheb_in[i] = what[i];
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
                    rhoz[n1][n2][i] = cheb_out[i]/2.0;
                }

                /* compute fourier coefficients of velocity derivative at mesh */
                /* points from uhat using FFT                                  */
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
            }
        }
    }
}

/******************************************************************************/
