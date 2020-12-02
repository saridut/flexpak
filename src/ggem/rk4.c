#include "rk4.h"
#include <stdlib.h>
#include <string.h>


//******************************************************************************

Rk4 rk4_create(const int nvar){

    Rk4    rk4;

    rk4.eos = true;
    rk4.info = 0;
    rk4.nvar = nvar;
    rk4.dt = 0.0;
    rk4.y0 = NULL;
    rk4.u = malloc(nvar*sizeof(double));
    rk4.f = malloc(nvar*sizeof(double));
    rk4.sol = malloc(nvar*sizeof(double));

    return rk4;
}

//******************************************************************************

void rk4_close(Rk4* this){

    this->y0 = NULL;

    free(this->u);
    free(this->f);
    free(this->sol);
}

//******************************************************************************

void rk4_set_initval(Rk4* this, double* y0){

    this->eos = false;
    this->info = 0;
    this->y0 = y0;
}

//******************************************************************************

void rk4_set_stepsize(Rk4* this, double dt){

    this->dt = dt;
}

//******************************************************************************

void rk4_integrate(Rk4* this){

    int     nvar = this->nvar;
    double  dt = this->dt;

    //Actions for this->info
    //0: first call to integrator
    //1-4: repeat call to integrator
    if (this->info == 0){
        for (int i=0; i<nvar; ++i){
            this->u[i] = this->y0[i];
            this->sol[i] = this->y0[i];
        }
        this->info = 1;
    }
    else if (this->info == 1){
        //f is f0. Calculating u1.
        for (int i=0; i<nvar; ++i){
            this->u[i] = this->y0[i] + this->f[i]*dt*0.5; //u1
            this->sol[i] +=  this->f[i]*(dt/6.0);
        }
        this->info = 2;
    }
    else if (this->info == 2){
        //f is f1. Calculating u2.
        for (int i=0; i<nvar; ++i){
            this->u[i] = this->y0[i] + this->f[i]*dt*0.5; //u2
            this->sol[i] +=  this->f[i]*(dt/3.0);
        }
        this->info = 3;
    }
    else if (this->info == 3){
        //f is f2. Calculating u3.
        for (int i=0; i<nvar; ++i){
            this->u[i] = this->y0[i] + this->f[i]*dt; //u3
            this->sol[i] +=  this->f[i]*(dt/3.0);
        }
        this->info = 4;
    }
    else if (this->info == 4){
        //f is f3. Calculating sol.
        for (int i=0; i<nvar; ++i){
            this->sol[i] +=  this->f[i]*(dt/6.0);
        }
        this->info = 0;
        this->eos = true;
    }

}

//******************************************************************************

void rk4_repeat(Rk4* this){

    /* Copy the solution to the initial value buffer.  Setting info=1 and   */
    /* eos=false to indicate integration should commence from computing u1. */
    /* Note that f is already set to f0 from previous step.                 */

    memcpy(this->y0, this->sol, (this->nvar)*sizeof(double));
    this->info = 1;
    this->eos = false;
}

//******************************************************************************

