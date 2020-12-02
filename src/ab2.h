#ifndef AB2_H
#define AB2_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

struct Ab2{
    int nvar;
    double  dt0; 
    double  dt1; 
    double* f0;
    double* y1;
    double* f1;
    double* sol;
};

typedef struct Ab2 Ab2;

Ab2 ab2_create(int nvar);

void ab2_close(Ab2* this);

void ab2_set_initval(Ab2* this, double* y1, double* f0);

void ab2_set_stepsize(Ab2* this, double dt1, double dt0);

void ab2_integrate(Ab2* this);

double* ab2_get_solution(Ab2* this);

void ab2_repeat(Ab2* this);

#ifdef __cplusplus
}
#endif

#endif // AB2_H
