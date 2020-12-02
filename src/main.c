#include "mesh.h"
#include "boundary_mesh.h"
#include "ref_config.h"
#include "mesh_transformation.h"
#include "ab2.h"
#include "assembler.h"
#include "free_space_solver.h"
#include "sim_params.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

int main()
{
    /*
     *  Read in simulation parameters
     */
    SimParams sp = sim_params_create();

    /*
     *  Create mesh
     */
    char fn_nodes[128];
    char fn_elements[128];
    char fn_conn[128];
    char fn_vtk[128];

    sprintf(fn_nodes, "%s/%s_nodes.txt", sp.model_dir, sp.model_name);
    sprintf(fn_elements, "%s/%s_elements.txt", sp.model_dir, sp.model_name);

    Mesh sheet = mesh_create(2, "tri");
    mesh_add_geometry_from_file(&sheet, fn_nodes);

    for (int i=0; i < 3; ++i){
        for (int j=0; j < 3; ++j){
            sprintf(fn_conn, "%s/%s_conn_%d-%d.txt", sp.model_dir,
                    sp.model_name, i, j);
            mesh_add_connectivity_from_file(&sheet, i, j, fn_conn);
        }
    }

    /*
     *  Create boundary mesh
     */
     BoundaryMesh bsheet = bmesh_create();
     bmesh_init(&bsheet, &sheet);
  
     /*
      *  Create reference configuration
      */
     RefConfig rsheet = refconf_create();
     refconf_init(&rsheet, &sheet);
  
     double dir[3] = {0.0, 0.0, 1.0};
     double axis[3] = {0.0, 1.0, 0.0};
     //mesh_rotate_axis_angle(&sheet, axis, M_PI_2);

     //mesh_perturb(&sheet, dir, 0.001);
     mesh_pull_back(&sheet);
     mesh_scale(&sheet, 1.5, 1.5, 1.5);
  
     /*
      *  Create assembler
      */
     Assembler as = asm_create();
     asm_init(&as, &sheet, &bsheet, &rsheet, sp.E, sp.kb);
  
     /*
      *  Create free-space solver
      */
     FreeSpaceSolver fss = fs_solver_create();
     fs_solver_add_mesh(&fss, &sheet);
     fs_solver_set_shear_rate(&fss, sp.gam_dot);
     fs_solver_set_num_quad_points(&fss, sp.num_quad_points[0],
             sp.num_quad_points[1]);
     fs_solver_set_sldq(&fss, as.traction);
     fs_solver_init(&fss);
    
     /*
      *  Create time integrator
      */
     int num_coordinates;
     double* velocity = NULL;
     double* coordinates = mesh_get_coordinates(&sheet, &num_coordinates);
     Ab2  timestepper = ab2_create(num_coordinates);
  
     ab2_set_initval(&timestepper, coordinates, velocity);
     ab2_set_stepsize(&timestepper, sp.dt, sp.dt);

    /* Time loop begins */
    double t = sp.t_start;

    int num_steps = 510;
    int istep = 0;

    /* Write output */
    sprintf(fn_vtk, "%s/sheet_t=%d.vtk", sp.output_dir, istep);
    mesh_to_vtk(&sheet, fn_vtk);

    while (istep < num_steps){
        printf("istep = %d\n", istep);
        asm_calc_sldq(&as);

        //printf("sldq~~~~~~~\n");
        //for(int i=0; i < 9; ++i){
        //    for (int j=0; j < 3; ++j){
        //        printf("%f  ", as.traction[3*i+j]);
        //    }
        //    printf("\n");
        //}

        fs_solver_solve(&fss);

        velocity = fs_solver_get_velocity(&fss);

        //printf("velocity~~~~~~~\n");
        //for(int i=0; i < 9; ++i){
        //    for (int j=0; j < 3; ++j){
        //        printf("%f  ", velocity[3*i+j]);
        //    }
        //    printf("\n");
        //}

        memcpy(timestepper.f1, velocity, num_coordinates*sizeof(double));

        ab2_integrate(&timestepper);

        double* coordinates = ab2_get_solution(&timestepper);

        mesh_set_coordinates(&sheet, coordinates);
        mesh_pull_back(&sheet);

        ab2_repeat(&timestepper);
        //t += sp.dt;
        istep+= 1;

        /* Write output */
        if (istep%1 == 0){
            sprintf(fn_vtk, "%s/sheet_t=%d.vtk", sp.output_dir, istep);
            mesh_to_vtk(&sheet, fn_vtk);
        }
    }


    /* Release resources */
    ab2_close(&timestepper);

    fs_solver_delete(&fss);

    refconf_delete(&rsheet);

    bmesh_delete(&bsheet);

    mesh_delete(&sheet);
}

