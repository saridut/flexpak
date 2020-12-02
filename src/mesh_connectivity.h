#ifndef MESH_CONNECTIVITY_H
#define MESH_CONNECTIVITY_H

#ifdef __cplusplus
extern "C" {
#endif

struct MeshConnectivity{
    /* Dimensions: Represents connections from entities of topological dimension */
    /* d0 to d1.                                                                 */
    int d0;
    int d1;
    int num_entities; /* Number of entities of topological dimension d0 */
    int num_connections; /* Total number of connections */
    int* connections;
    int* ind_pos;    /* Index array */
};

typedef struct MeshConnectivity MeshConnectivity;

/* Creates an empty MeshConnectivity from entities of topological dimension
 * d0 to entities of topological dimension d1. */
MeshConnectivity meshcon_create(const int d0, const int d1);


/* Initializes from array */
void meshcon_from_array(MeshConnectivity* this, const int num_entities,
        const int num_connections, const int* connections, const int* ind_pos);


/* Initializes by reading connectivity data from file */
void meshcon_from_file(MeshConnectivity* this, const char* fn);


/* Copies a MeshConnectivity object of the same topological dimension pair */
/* (d0, d1)                                                                */
void meshcon_copy(MeshConnectivity* this, const MeshConnectivity* other);


/* Delete a MeshConnectivity object */
void meshcon_delete(MeshConnectivity* this);


/* Clears an existing MeshConnectivity object */
void meshcon_clear(MeshConnectivity* this);


/* Returns the number of entities of topological dimension d0 */
int meshcon_get_num_entities(const MeshConnectivity* this);


/* Returns the total number of connections in (d0->d1) */
int meshcon_get_num_connections(const MeshConnectivity* this);


/* Returns the number of connections for entity ient */
int meshcon_get_num_conns(const MeshConnectivity* this, const int ient);


/* Gets a pointer to connections for entity ient */
int* meshcon_get_conns(const MeshConnectivity* this, const int ient,
        int* num_conns);


/* Gets a pointer to connection i for entity ient */
int* meshcon_get_conn(const MeshConnectivity* this, const int ient,
        const int i);


/* Writes to a text file */
void meshcon_to_file(MeshConnectivity* this, const char* fn);

#ifdef __cplusplus
}
#endif

#endif //MESH_CONNECTIVITY_H


