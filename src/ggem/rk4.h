#ifndef RK4_H
#define RK4_H

#include <stdbool.h>

struct Rk4{
    bool eos;
    int info;
    int nvar;
    double  dt;
    double* y0;
    double* u;
    double* f;
    double* sol;
};

typedef struct Rk4 Rk4;

Rk4 rk4_create(int nvar);

void rk4_close(Rk4* this);

void rk4_set_initval(Rk4* this, double* y0);

void rk4_set_stepsize(Rk4* this, double dt);

void rk4_integrate(Rk4* this);

void rk4_repeat(Rk4* this);

#endif // RK4_H
