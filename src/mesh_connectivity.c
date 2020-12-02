#include "mesh_connectivity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/******************************************************************************/

/* Creates an empty MeshConnectivity from entities of topological dimension
 * d0 to entities of topological dimension d1. */
MeshConnectivity meshcon_create(const int d0, const int d1){

    MeshConnectivity mc;

    mc.d0 = d0;
    mc.d1 = d1;
    mc.num_entities = 0;
    mc.num_connections = 0;
    mc.connections = NULL;
    mc.ind_pos = NULL;

    return mc;
}

/******************************************************************************/

/* Initialize from array */
void meshcon_from_array(MeshConnectivity* this, const int num_entities,
        const int num_connections, const int* connections, const int* ind_pos){

    this->num_entities = num_entities;
    this->num_connections = num_connections;

    this->ind_pos = malloc((num_entities+1)*sizeof(int));
    this->connections = malloc(num_connections*sizeof(int));

    memcpy(this->ind_pos, ind_pos, (num_entities+1)*sizeof(int));
    memcpy(this->connections, connections, num_connections*sizeof(int));
}

/******************************************************************************/

/* Initializes by reading connectivity data from file */
void meshcon_from_file(MeshConnectivity* this, const char* fn){

    FILE*    fp = fopen(fn, "r");
    char*    line=NULL; /* pointer to a line buffer */
    size_t   n=0;       /* size of line buffer */
    char*    word;
    int      nc;        /* number of connections for entity ient */
    int      beg;

    /* Read the header line */
    getline(&line, &n, fp);
    this->num_entities = strtol(strtok(line, " "), NULL, 10);
    this->num_connections = strtol(strtok(NULL, " "), NULL, 10);

    /* Allocate memory for arrays ind_pos and connections */
    this->ind_pos = malloc((this->num_entities+1)*sizeof(int));
    this->connections = malloc(this->num_connections*sizeof(int));

    /* Loop over all lines. For each line read a sequence of ints -- the */
    /* first int corresponds to the index of the connection.             */
    this->ind_pos[0] = 0;
    for (int ient=0; ient < this->num_entities; ++ient){
        nc = 0; /* number of connections for ient */
        beg = this->ind_pos[ient];
        getline(&line, &n, fp);
        word = strtok(line, " "); /* Discard the connection index */
        while (true) {
            word = strtok(NULL, " ");
            if (word != NULL){
                this->connections[beg+nc] = strtol(word, NULL, 10);
                nc += 1;
            }
            else {
                break;
            }
        }
        this->ind_pos[ient+1] = this->ind_pos[ient] + nc;
    }

    fclose(fp);
}

/******************************************************************************/

/* Copies a MeshConnectivity object of the same topological dimension pair */
/* (d0, d1)                                                                */
void meshcon_copy(MeshConnectivity* this, const MeshConnectivity* other){

    int tnc = other->num_connections;
    int ne = other->num_entities;

    /* Allocate memory for new arrays */
    this->connections = realloc(this->connections, tnc*sizeof(int));
    this->ind_pos = realloc(this->ind_pos, (ne+1)*sizeof(int));

    /* Copy from other */
    this->num_connections = tnc;
    this->num_entities = ne;

    memcpy(this->connections, other->connections, tnc*sizeof(int));
    memcpy(this->ind_pos, other->ind_pos, (ne+1)*sizeof(int));
}

/******************************************************************************/

/* Delete a MeshConnectivity object */
void meshcon_delete(MeshConnectivity* this){

    free(this->connections) ;
    free(this->ind_pos) ;

    this->num_entities = 0;
    this->num_connections = 0;
}

/******************************************************************************/

/* Clears an existing MeshConnectivity object */
void meshcon_clear(MeshConnectivity* this){

    free(this->connections);
    free(this->ind_pos);

    this->ind_pos = NULL;
    this->connections = NULL;

    this->num_entities = 0;
    this->num_connections = 0;
}

/******************************************************************************/

/* Returns the number of entities of topological dimension d0 */
int meshcon_get_num_entities(const MeshConnectivity* this){

    return this->num_entities;
}

/******************************************************************************/

/* Returns the total number of connections in (d0->d1) */
int meshcon_get_num_connections(const MeshConnectivity* this){

    return this->num_connections;
}

/******************************************************************************/

/* Returns the number of connections for entity ient */
int meshcon_get_num_conns(const MeshConnectivity* this, const int ient){

    return (this->ind_pos[ient+1] - this->ind_pos[ient]);
}

/******************************************************************************/

/* Gets a pointer to connections for entity ient */
int* meshcon_get_conns(const MeshConnectivity* this, const int ient,
        int* num_conns){

    int beg = this->ind_pos[ient];

    if (num_conns != NULL){
        *num_conns = this->ind_pos[ient+1] - this->ind_pos[ient];
    }
    return  (this->connections + beg);
}

/******************************************************************************/

/* Gets a pointer to connection i for entity ient */
int* meshcon_get_conn(const MeshConnectivity* this, const int ient,
        const int i){

    int beg = this->ind_pos[ient];

    return (this->connections + beg + i);
}

/******************************************************************************/

/* Writes to a text file */
void meshcon_to_file(MeshConnectivity* this, const char* fn){

    FILE* fp = fopen(fn, "w");
    int   beg;
    int   end;

    /* Writing the header */
    fprintf(fp, "%d  %d\n", this->num_entities, this->num_connections);

    for (int ient=0; ient  < this->num_entities; ++ient){
        fprintf(fp, "%d", ient);

        beg = this->ind_pos[ient];
        end = this->ind_pos[ient+1];
        
        for (int j=beg; j<end; ++j){
            fprintf(fp, "    %d", this->connections[j]);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
}

/******************************************************************************/
