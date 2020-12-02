#include "mesh_geometry.h"
#include "mesh_connectivity.h"
#include <stdio.h>

int main(){

    const char* fn = "conn_2-0.txt";

//    MeshGeometry mg = meshgeom_create();
//
//    meshgeom_from_file(&mg, fn);
//
//    printf("Number of vertices %d\n", meshgeom_get_num_vertices(&mg));
//    printf("Number of coordinates %d\n", meshgeom_get_num_coordinates(&mg));
//
//    meshgeom_to_file(&mg, "to_file.txt");
//
//    meshgeom_delete(&mg);

    const int* conns=NULL;
    int nc;

    MeshConnectivity mc = meshcon_create(2, 0);

    meshcon_from_file(&mc, fn);

    printf("Number of entities %d\n", meshcon_get_num_entities(&mc));

    nc = meshcon_get_connections(&mc, 0, &conns);

    printf("Number of connections %d\n", nc);
    for (int i=0; i<nc; ++i){
        printf("%d  ", conns[i]);
    }

    printf("\n");

    meshcon_to_file(&mc, "to_file.txt");

    meshcon_delete(&mc);


}
