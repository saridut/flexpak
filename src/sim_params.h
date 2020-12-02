#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

#include <stdbool.h>

struct SimParams{
    char  model_dir[80];
    char  model_name[80];
    char  output_dir[80];

    double  t_start;
    double  t_end;
    double  dt;

    double  gam_dot;
    int  num_quad_points[2];
    double  ca; //capillary number
    double  E;  //initial Young's modulus
    double  kb; //bending mmodulus
};

typedef struct SimParams SimParams;

SimParams sim_params_create();

void sim_params_delete(SimParams* this);

void sim_params_from_file(SimParams* this);

void sim_params_to_file(SimParams* this, const char* fn);

void sim_params_check();

bool sim_params_is_comment(const char* linebuf);

void sim_params_read_value(SimParams* this, char* linebuf);

#endif //SIM_PARAMS_H
