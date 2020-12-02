#include "ab2.h"
#include <utils_math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//******************************************************************************

Ab2 ab2_create(const int nvar){

    Ab2    ab2;

    ab2.nvar = nvar;
    ab2.dt0 = 0.0;
    ab2.dt1 = 0.0;
    ab2.f0 = malloc(nvar*sizeof(double));
    ab2.y1 = malloc(nvar*sizeof(double));
    ab2.f1 = malloc(nvar*sizeof(double));
    ab2.sol = malloc(nvar*sizeof(double));

    return ab2;
}

//******************************************************************************

void ab2_close(Ab2* this){

    free(this->f0);
    free(this->y1);
    free(this->f1);
    free(this->sol);
}

//******************************************************************************

void ab2_set_initval(Ab2* this, double* y1, double* f0){

    if (f0 == NULL){
        zero_out(this->nvar, this->f0);
    }
    else {
        memcpy(this->f0, f0, this->nvar*sizeof(double));
    }
    memcpy(this->y1, y1, this->nvar*sizeof(double));
}

//******************************************************************************

void ab2_set_stepsize(Ab2* this, double dt1, double dt0){

    this->dt0 = dt0;
    this->dt1 = dt1;
}

//******************************************************************************

void ab2_integrate(Ab2* this){

    int     nvar = this->nvar;
    double  dt0 = this->dt0; 
    double  dt1 = this->dt1; 

    //Evaluate f(y1) and put it in f1
    for (int i=0; i<nvar; ++i){
        this->sol[i] = this->y1[i] + (0.5*dt1)*(3.0*this->f1[i] - this->f0[i]);
    }
}

//*****************************************************************************

double* ab2_get_solution(Ab2* this){

    return this->sol;
}

//*****************************************************************************

void ab2_repeat(Ab2* this){
    
    memcpy(this->f0, this->f1, (this->nvar)*sizeof(double));
    memcpy(this->y1, this->sol, (this->nvar)*sizeof(double));
}

//*****************************************************************************
