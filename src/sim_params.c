#include "sim_params.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*****************************************************************************/

SimParams sim_params_create()
{
    SimParams sp;

    sim_params_from_file(&sp);

    return sp;
}

/****************************************************************************/

void sim_params_delete(SimParams* this)
{
    /*
     *  Not implemented
     */
    ;
}

/****************************************************************************/

void sim_params_from_file(SimParams* this)
{
    FILE* fp = fopen("params.txt", "r");
    char* linebuf = NULL;
    size_t buflen = 0;

    while (getline(&linebuf, &buflen, fp) > 0 ){
        //Handling newline at the end
        linebuf[strcspn(linebuf, "\n")] = '\0';
        if (sim_params_is_comment(linebuf)){
            continue;
        }
        else {
            sim_params_read_value(this, linebuf);
        }
    }

    if (linebuf != NULL){
        free(linebuf); /* allocated initially by getline */
    }
    fclose(fp);
}

/****************************************************************************/

void sim_params_to_file(SimParams* this, const char* fn)
{
    /*
     *  Not Implemented
     */
    ;
}

/****************************************************************************/

void sim_params_check()
{
    /*
     *  Not implemented
     */
    ;
}

/*****************************************************************************/

/* A line is a comment if it is either blank or its leading non-blank
 * character is ``#''.
 */
bool sim_params_is_comment(const char* linebuf)
{
    // const char comchar[3] = {'#', ' ', '\t'};
    bool ret = true;

    while (*linebuf != '\0'){
        char c = *linebuf;
        if ( (c != '#') && (c != ' ') && (c != '\t') ) {
            ret = false;
            break;
        }
        ++linebuf;
    }
    return ret;
}

/*****************************************************************************/

void sim_params_read_value(SimParams* this, char* linebuf)
{
    char* key = strtok(linebuf, " \t");
    char* value = strtok(NULL, " \t");

    if (strcmp(key, "model_dir") == 0) {
        strcpy(this->model_dir, value);
    }
    else if (strcmp(key, "model_name") == 0) {
        strcpy(this->model_name, value);
    }
    else if (strcmp(key, "output_dir") == 0) {
        strcpy(this->output_dir, value);
    }
    else if (strcmp(key, "t_start") == 0) {
        this->t_start = strtod(value, NULL);
    }
    else if (strcmp(key, "t_end") == 0) {
        this->t_end = strtod(value, NULL);
    }
    else if (strcmp(key, "dt") == 0) {
        this->dt = strtod(value, NULL);
    }
    else if (strcmp(key, "gam_dot") == 0) {
        this->gam_dot = strtod(value, NULL);
    }
    else if (strcmp(key, "E") == 0) {
        this->E = strtod(value, NULL);
    }
    else if (strcmp(key, "kb") == 0) {
        this->kb = strtod(value, NULL);
    }
    else if (strcmp(key, "num_quad_points") == 0) {
        this->num_quad_points[0] = strtol(value, NULL, 10);
        char* value = strtok(NULL, " \t");
        this->num_quad_points[1] = strtol(value, NULL, 10);
    }
    else{
        printf("Unknown value %s\n", key);
    }
}

/*****************************************************************************/
