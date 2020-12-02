#include "assembler.h"
#include "triangle_cell.h"
#include "rotation.h"
#include "utils_math.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/******************************************************************************/

Assembler asm_create()
{
    Assembler asmb;

    asmb.msh = NULL;
    asmb.bmsh = NULL;
    asmb.refconf = NULL;
    asmb.cell_center = NULL;
    asmb.cell_area = NULL;
    asmb.cell_angle = NULL;
    asmb.cell_angle_cot = NULL;
    asmb.cell_normal = NULL;
    asmb.edge_len_sq = NULL;
    asmb.edge_center = NULL;
    asmb.amixed = NULL;
    asmb.vert_normal = NULL;
    asmb.mean_curv = NULL;
    asmb.gauss_curv = NULL;
    asmb.sl_op = NULL;
    asmb.sgmc = NULL;
    asmb.slmc = NULL;
    asmb.is_obtuse = NULL;
    asmb.traction = NULL;
    asmb.E = 0.0;
    asmb.kb = -1.0;

    return asmb;
}

/******************************************************************************/

/* Allocate memory and setup buffers for force assembly. If bending */
/* effects are not required pass a negative value for kb. If there  */
/* is no boundary, pass a null pointer as bmsh.                     */
void asm_init(Assembler* this, Mesh* msh, BoundaryMesh* bmsh,
        RefConfig* refconf, double E, double kb)
{
    this->msh = msh;
    this->refconf = refconf;
    this->bmsh = bmsh;
    this->E = E;
    this->kb = kb;

    int num_cells = mesh_get_num_entities(this->msh, 2);
    int num_edges = mesh_get_num_entities(this->msh, 1);
    int num_vertices = mesh_get_num_entities(this->msh, 0);

    this->cell_center = malloc(3*num_cells*sizeof(double));
    this->cell_angle = malloc(3*num_cells*sizeof(double));
    this->cell_angle_cot = malloc(3*num_cells*sizeof(double));
    this->cell_normal = malloc(3*num_cells*sizeof(double));
    this->cell_area = malloc(num_cells*sizeof(double));
    this->is_obtuse = malloc(num_cells*sizeof(int));

    this->edge_len_sq = malloc(num_edges*sizeof(double));
    this->edge_center = malloc(3*num_edges*sizeof(double));

    this->amixed = malloc(num_vertices*sizeof(double));
    this->mean_curv = malloc(num_vertices*sizeof(double));
    this->gauss_curv = malloc(num_vertices*sizeof(double));
    this->vert_normal = malloc(3*num_vertices*sizeof(double));

    this->sl_op = malloc(3*num_vertices*sizeof(double));
    this->sgmc = malloc(3*num_vertices*sizeof(double));
    this->slmc = malloc(num_vertices*sizeof(double));
    this->traction = malloc(3*num_vertices*sizeof(double));

}

/******************************************************************************/

void asm_delete(Assembler* this)
{

    this->msh = NULL;
    this->bmsh = NULL;
    this->refconf = NULL;
    free (this->cell_center);
    free (this->cell_area);
    free (this->cell_angle);
    free (this->cell_angle_cot);
    free (this->cell_normal);
    free (this->edge_len_sq);
    free (this->edge_center);
    free (this->amixed);
    free (this->vert_normal);
    free (this->mean_curv);
    free (this->gauss_curv);
    free (this->sl_op);
    free (this->sgmc);
    free (this->slmc);
    free (this->is_obtuse);
    free(this->traction);
}

/******************************************************************************/

/*
 *  TODO: As of now, will overwrite traction. MUST FIX THIS.
 */
void asm_calc_sldq(Assembler* this)
{
    int num_coordinates = mesh_get_num_coordinates(this->msh);
    asm_calc_traction(this);
    for (int i=0; i < num_coordinates; ++i){
        this->traction[i] /= -(4.0*M_PI);
    }
}

/******************************************************************************/

void asm_calc_traction(Assembler* this)
{
    asm_calc_traction_tension(this);

    if (this->kb > 0){
        asm_calc_traction_bending(this);
    }

    if (this->bmsh != NULL){
        asm_set_bc(this);
    }
}

/******************************************************************************/

/* Returns a pointer to the traction vector at all nodes */
double* asm_get_traction(const Assembler* this, int* n)
{
    if (n != NULL){
        int num_vertices = mesh_get_num_entities(this->msh, 0);
        *n = 3*num_vertices;
    }
    return this->traction;
}

/******************************************************************************/

/* Returns a pointer to the traction vector at vert ivert */
double* asm_get_traction_at_vert(const Assembler* this, const int ivert)
{
    return (this->traction+3*ivert);
}

/******************************************************************************/

/* Returns a pointer to the single layer density at all nodes
 * TODO: As of now, this returns the overwritten traction. FIX THIS */
double* asm_get_sldq(const Assembler* this, int* n)
{
    if (n != NULL){
        int num_vertices = mesh_get_num_entities(this->msh, 0);
        *n = 3*num_vertices;
    }
    return this->traction;
}

/******************************************************************************/

void asm_calc_traction_tension(Assembler* this)
{
    double u[3];
    double v[3];

    double dG11_du[3];
    double dG11_dv[3];

    double dG22_du[3];
    double dG22_dv[3];

    double dG12_du[3];
    double dG12_dv[3];

    double dlambda1_du[3];
    double dlambda1_dv[3];

    double dlambda2_du[3];
    double dlambda2_dv[3];

    double g_cell_coords[9];
    double l_cell_coords[9];
    double dcm[9];
    double Aop[9];
    double Bop[9];
    double C[9];
    double g_Fp[9];
    double l_Fp[9];

    double* A;
    double* B;
    double* coords;
    int*    conns;

    double cell_area;
    double Ac;
    double G11, G22, G12;
    double du_dx, du_dy, dv_dx, dv_dy;
    double du_dxp1, dv_dyp1;
    double tmp_a, tmp_b;
    double lambda1, lambda2;
    double dW_dlambda1, dW_dlambda2;
    double rlambda1, rlambda2, rlambda12, rlsq;

    int nc;
    int num_cells = mesh_get_num_entities(this->msh, 2);
    int num_verts = mesh_get_num_vertices(this->msh);

    /* Zeroing out traction buffer */
    zero_out(3*num_verts, this->traction);


    //printf("traction before traction\n");
    //for(int i=0; i < num_verts; ++i){
    //    for (int j=0; j < 3; ++j){
    //        printf("%f  ", this->traction[3*i+j]);
    //    }
    //    printf("\n");
    //}
    for (int icell=0; icell < num_cells; ++icell){
        mesh_get_ent_coords(this->msh, 2, icell, g_cell_coords, NULL);

        /* Calculating area of all cells (i.e. triangles) */
        this->cell_area[icell] = tri_area(g_cell_coords);

        /* Converting vertices to local coordinates; getting the dcm as well */
        tri_to_local(g_cell_coords, l_cell_coords, dcm);

        //printf("l_cell_coords\n");
        //for(int i=0; i < 3; ++i){
        //    for (int j=0; j < 3; ++j){
        //        printf("%f  ", l_cell_coords[3*i+j]);
        //    }
        //    printf("\n");
        //}

        /* Calculating displacement of the vertices */
        zero_out(3, u);
        zero_out(3, v);

        coords = refconf_get_coords(this->refconf, icell, NULL);

        u[1] = l_cell_coords[3*1+0] - coords[0];
        u[2] = l_cell_coords[3*2+0] - coords[1];
        v[2] = l_cell_coords[3*2+1] - coords[2];

        //printf("u: %f  %f  %f\n", u[0], u[1], u[2]);
        //printf("v: %f  %f  %f\n", v[0], v[1], v[2]);

        /* Calculating derivatives of G w.r.t. u and v */
        A = refconf_get_A(this->refconf, icell, NULL);
        B = refconf_get_B(this->refconf, icell, NULL);

        outer(3, 3, A, A, Aop);
        outer(3, 3, B, B, Bop);

        /* C = B*A^T + A*B^T */
        mm_mult(3, 1, B, 3, A, C);
        add_transpose(3, C);

        du_dx = dot(3, u, A);
        du_dy = dot(3, u, B);
        dv_dx = dot(3, v, A);
        dv_dy = dot(3, v, B);

        du_dxp1 = du_dx + 1.0;
        dv_dyp1 = dv_dy + 1.0;

        G11 = du_dxp1*du_dxp1 + dv_dx*dv_dx;
        G22 = du_dy*du_dy + dv_dyp1*dv_dyp1;
        G12 = du_dxp1*du_dy + dv_dyp1*dv_dx;

        tmp_a = G11 + G22;
        tmp_b = sqrt((G11-G22)*(G11-G22) + 4.0*G12*G12);
        lambda1 = sqrt(0.5*(tmp_a + tmp_b));
        lambda2 = sqrt(0.5*(tmp_a - tmp_b));

        mv_mult(3, 3, Aop, u, dG11_du);
        mv_mult(3, 3, Aop, v, dG11_dv);
        mv_mult(3, 3, Bop, u, dG22_du);
        mv_mult(3, 3, Bop, v, dG22_dv);
        mv_mult(3, 3, C, u, dG12_du);
        mv_mult(3, 3, C, v, dG12_dv);

        for (int i=0; i < 3; ++i){
            dG11_du[i] += A[i];
            dG22_dv[i] += B[i];
            dG12_du[i] += B[i];
            dG12_dv[i] += A[i];
        }

        for (int i=0; i < 3; ++i){
            dG11_du[i] *= 2;
            dG11_dv[i] *= 2;
            dG22_du[i] *= 2;
            dG22_dv[i] *= 2;
        }

        rlambda1 = 1.0/lambda1;
        rlambda2 = 1.0/lambda2;
        rlambda12 = rlambda1*rlambda2;
        rlsq = rlambda12*rlambda12;
        dW_dlambda1 = this->E*(lambda1 - rlambda1*rlsq)/3.0;
        dW_dlambda2 = this->E*(lambda2 - rlambda2*rlsq)/3.0;
        //printf("lambda1: %f\n", lambda1);
        //printf("lambda2: %f\n", lambda2);
        //printf("dW_dlambda1: %f\n", dW_dlambda1);
        //printf("dW_dlambda2: %f\n", dW_dlambda2);
        //printf("tmp_b: %f\n", tmp_b);
        //printf("dG11_du: %f  %f  %f\n", dG11_du[0], dG11_du[1], dG11_du[2]);
        //printf("dG22_du: %f  %f  %f\n", dG22_du[0], dG22_du[1], dG22_du[2]);

        if (isclose(tmp_b, 0.0, 1e-12, 0.0)){
            for (int i=0; i < 3; ++i){
                dlambda1_du[i] = (dG11_du[i] + dG22_du[i])/(4.0*lambda1);
                dlambda1_dv[i] = (dG11_dv[i] + dG22_dv[i])/(4.0*lambda1);
                dlambda2_du[i] = (dG11_du[i] + dG22_du[i])/(4.0*lambda2);
                dlambda2_dv[i] = (dG11_dv[i] + dG22_dv[i])/(4.0*lambda2);
            }
        }
        else {
            for (int i=0; i < 3; ++i){
                dlambda1_du[i] = (dG11_du[i] + dG22_du[i] 
                        + ((G11-G22)*(dG11_du[i]-dG22_du[i])
                        + 4.0*G12*dG12_du[i])/tmp_b)/(4.0*lambda1);

                dlambda1_dv[i] = (dG11_dv[i] + dG22_dv[i] 
                        + ((G11-G22)*(dG11_dv[i]-dG22_dv[i])
                        + 4.0*G12*dG12_dv[i])/tmp_b)/(4.0*lambda1);

                dlambda2_du[i] = (dG11_du[i] + dG22_du[i] 
                        - ((G11-G22)*(dG11_du[i]-dG22_du[i])
                        + 4.0*G12*dG12_du[i])/tmp_b)/(4.0*lambda2);

                dlambda2_dv[i] = (dG11_dv[i] + dG22_dv[i] 
                        - ((G11-G22)*(dG11_dv[i]-dG22_dv[i])
                        + 4.0*G12*dG12_dv[i])/tmp_b)/(4.0*lambda2);
            }
        }

        cell_area = this->cell_area[icell];
        for (int i=0; i < 3; ++i){
            l_Fp[3*i+0] = cell_area*(dW_dlambda1*dlambda1_du[i]
                        + dW_dlambda2*dlambda2_du[i]);
            l_Fp[3*i+1] = cell_area*(dW_dlambda1*dlambda1_dv[i]
                        + dW_dlambda2*dlambda2_dv[i]);
            l_Fp[3*i+2] = 0.0;
        }

        //for (int k=0; k < 3; ++k){
        //    printf("l_Fp: %f  %f  %f\n", l_Fp[3*k+0], l_Fp[3*k+1], l_Fp[3*k+2]);
        //}
        /* Shift to global basis */
        shift_vectors_dcm(3, l_Fp, dcm, true, g_Fp);

        //for (int k=0; k < 3; ++k){
        //    printf("g_Fp: %f  %f  %f\n", g_Fp[3*k+0], g_Fp[3*k+1], g_Fp[3*k+2]);
        //}
        /* Add reaction forces at each vertex incident to icell */
        conns = mesh_get_conns(this->msh, 2, 0, icell, &nc);
        for (int i=0; i< nc; ++i){
            for (int j=0; j < 3; ++j){
                this->traction[3*conns[i]+j] += g_Fp[3*i+j];
            }
        }
    }

    /* Area weighted average of the reaction forces */
    for (int ivert=0; ivert < num_verts; ++ivert){
        conns = mesh_get_conns(this->msh, 0, 2, ivert, &nc);
        Ac = 0.0;
        for (int i=0; i< nc; ++i){
            Ac += this->cell_area[conns[i]];
        }
        for (int i=0; i < 3; ++i){
            /* No negative sign here. Is Charrier et al's force directed */
            /* in opposite direction?                                    */
            this->traction[3*ivert+i] /= (Ac/3.0);
        }
    }

    //printf("traction (tension) =======\n");
    //for(int i=0; i < 9; ++i){
    //    for (int j=0; j < 3; ++j){
    //        printf("%f  ", this->traction[3*i+j]);
    //    }
    //    printf("\n");
    //}

}

/******************************************************************************/

void asm_calc_traction_bending(Assembler* this)
{
    int num_cells = mesh_get_num_entities(this->msh, 2);
    double g_cell_coords[9];

    /* Calculating normals of all cells (i.e. triangles) */
    for (int icell=0; icell < num_cells; ++icell){
        mesh_get_ent_coords(this->msh, 2, icell, g_cell_coords, NULL);
        tri_normal(g_cell_coords, this->cell_normal+3*icell);
    }

    /* Calculate the length of each edge */
    asm_calc_edge_len_sq(this);

    /* Calculating angles of each cell and their cotangents */
    asm_calc_angles(this);

    //printf("angles @@@@@@@@@@\n");
    //for (int l=0; l < num_cells; ++l){
    //    for (int m=0; m < 3; ++m){
    //        printf("%f ", this->cell_angle[3*l+m]);
    //    }
    //    printf("\n");
    //}

    /* Calculating edge centers */
    asm_calc_edge_centers(this);

    /* Calculating center of cells -- circumcenter for non-obtuse */
    /* triangles, else the midpoint of the side opposite to the   */
    /* obtuse angle                                               */
    asm_calc_cell_centers(this);

    /* Calculating the mixed area for each vertex */
    int num_verts = mesh_get_num_vertices(this->msh);

    for (int ivert=0; ivert < num_verts; ++ivert){
        double amixed = 0.0;
        double gauss_curv = 2*M_PI;
        /* conns_tri is a pointer to triangles in the 1-ring neighborhood */
        /* of vertex ivert                                                */
        int nc_tri;
        int l_e0, l_e1;
        int g_e0, g_e1;
        int g_v0, g_v1;
        int* conns_tri = mesh_get_conns(this->msh, 0, 2, ivert, &nc_tri);
        double  *v0_coords, *v1_coords, *v_coords;
        double sl_op[3] = {0.0, 0.0, 0.0};

        for (int j=0; j < nc_tri; ++j){
            int icell = conns_tri[j];
            int l_ivert = mesh_get_local_index(this->msh, 2, icell, 0, ivert);
            switch (l_ivert){
                case 0:
                    l_e0 = 1;
                    l_e1 = 2;
                    break;
                case 1:
                    l_e0 = 0;
                    l_e1 = 2;
                    break;
                case 2:
                    l_e0 = 0;
                    l_e1 = 1;
                    break;
            }
            g_e0 = mesh_get_global_index(this->msh, 2, icell, 1, l_e0);
            g_e1 = mesh_get_global_index(this->msh, 2, icell, 1, l_e1);

            g_v0 = mesh_get_global_index(this->msh, 2, icell, 0, l_e0);
            g_v1 = mesh_get_global_index(this->msh, 2, icell, 0, l_e1);

            v0_coords = mesh_get_coords(this->msh, g_v0, NULL);
            v1_coords = mesh_get_coords(this->msh, g_v1, NULL);
            v_coords = mesh_get_coords(this->msh, ivert, NULL);

            for (int i=0; i < 3; ++i){
                sl_op[i] += (v0_coords[i]-v_coords[i])*
                            this->cell_angle_cot[3*icell+l_e0] +
                            (v1_coords[i]-v_coords[i])*
                            this->cell_angle_cot[3*icell+l_e1];
            }

            if (this->is_obtuse[icell] == -1){
                /* Apply Voronoi formula */
                amixed += (this->edge_len_sq[g_e0]*
                            this->cell_angle_cot[3*icell+l_e0] +
                            this->edge_len_sq[g_e1]*
                            this->cell_angle_cot[3*icell+l_e1])/8.0;
            }
            else {
                /* angle at ivert is obtuse */
                if (this->is_obtuse[icell] == l_ivert){
                    amixed += this->cell_area[icell]/2.0;
                }
                /* angle at ivert is not obtuse */
                else {
                    amixed += this->cell_area[icell]/4.0;
                }
            }

            gauss_curv -= this->cell_angle[3*icell+l_ivert];
            //if (ivert==4) printf("%f  %f\n", this->cell_angle[3*icell+l_ivert], gauss_curv);
        }
        memcpy(this->sl_op+3*ivert, sl_op, 3*sizeof(double));
        this->amixed[ivert] = amixed;
        this->gauss_curv[ivert] = gauss_curv;
    }

    for (int ivert=0; ivert < num_verts; ++ivert){
        double amixed = this->amixed[ivert];
        for (int i=0; i < 3; ++i){
            this->sl_op[3*ivert+i] /= (2.0*amixed);
        }
        double norm_sl_op = norm(3, this->sl_op+3*ivert);
        this->mean_curv[ivert] = norm_sl_op/2.0;
        for (int i=0; i < 3; ++i){
            int m = 3*ivert + i;
            this->vert_normal[m] = this->sl_op[m]/norm_sl_op;
        }
        this->gauss_curv[ivert] /= amixed;
    }

    /* Calculation for the surface gradient of mean curvature */

    asm_calc_sgmc(this);

    /* Integration of surface divergence over the contour of vertex patch */
    for (int ivert=0; ivert < num_verts; ++ivert){
        int nc_tri;
        double* vert = mesh_get_coords(this->msh, ivert, NULL);
        int* conns_tri = mesh_get_conns(this->msh, 0, 2, ivert, &nc_tri);
        double* cell_center;
        double* cell_normal;
        int l_ivert;
        int l_e0, l_e1;
        int g_e0, g_e1;
        double *ec_0, *ec_1;
        double p0[3], p1[3], r[3];
        double p0_hat[3], q0[3];
        double p1_hat[3], q1[3];
        double len_p0, len_p1;

        for (int j=0; j < nc_tri; ++j){
            int icell = conns_tri[j];
            /* pointer to cell center coordinates */
            cell_center = this->cell_center+3*icell;
            cell_normal = this->cell_normal+3*icell;
            l_ivert = mesh_get_local_index(this->msh, 2, icell, 0, ivert);
            switch (l_ivert){
                case 0:
                    l_e0 = 1;
                    l_e1 = 2;
                    break;
                case 1:
                    l_e0 = 0;
                    l_e1 = 2;
                    break;
                case 2:
                    l_e0 = 0;
                    l_e1 = 1;
                    break;
            }
            g_e0 = mesh_get_global_index(this->msh, 2, icell, 1, l_e0);
            g_e1 = mesh_get_global_index(this->msh, 2, icell, 1, l_e1);
            ec_0 = this->edge_center+3*g_e0;
            ec_1 = this->edge_center+3*g_e1;
            for (int i=0; i < 3; ++i){
                p0[i] = ec_0[i] - cell_center[i];
                p1[i] = ec_1[i] - cell_center[i];
                r[i]  = vert[i] - cell_center[i];
            }
            len_p0 = norm(3, p0);
            len_p1 = norm(3, p1);

            if (!isclose(len_p0, 0.0, 1e-8, 1e-10)){
                for (int i=0; i < 3; ++i){
                    p0_hat[i] = p0[i]/len_p0;
                }
                cross(cell_normal, p0_hat, q0);
                unitize(3, q0);
                if (dot(3, q0, r) > 0.0){
                    q0[0] = -q0[0]; q0[1] = -q0[1]; q0[2] = -q0[2];
                }
                this->slmc[ivert] += dot(3, &this->sgmc[icell], q0)*len_p0;
            }

            if (!isclose(len_p1, 0.0, 1e-8, 1e-10)){
                for (int i=0; i < 3; ++i){
                    p1_hat[i] = p1[i]/len_p1;
                }
                cross(cell_normal, p1_hat, q1);
                unitize(3, q1);
                if (dot(3, q1, r) > 0.0){
                    q1[0] = -q1[0]; q1[1] = -q1[1]; q1[2] = -q1[2];
                }
                this->slmc[ivert] += dot(3, &this->sgmc[icell], q1)*len_p1;
            }
        }
        this->slmc[ivert] /= this->amixed[ivert];
    }

    /* Calculation for traction due to bending */
    double mean_curv;
    double gauss_curv;
    double slmc;
    for (int ivert=0; ivert < num_verts; ++ivert){
        mean_curv = this->mean_curv[ivert];
        gauss_curv = this->gauss_curv[ivert];
        slmc = this->slmc[ivert];
        for (int i=0; i < 3; ++i){
            this->traction[3*ivert+i] += 
                 -this->kb*(2*slmc + 4*mean_curv*(mean_curv*mean_curv 
                    - gauss_curv))*this->vert_normal[3*ivert+i];
        }
    }
    //printf("traction (bending) ~~~~~~~\n");
    //for (int ivert=0; ivert < num_verts; ++ivert){
    //    mean_curv = this->mean_curv[ivert];
    //    gauss_curv = this->gauss_curv[ivert];
    //    slmc = this->slmc[ivert];
    //    //if (ivert == 4){
    //    //    printf("%f %f %f ", this->amixed[ivert], gauss_curv, mean_curv);
    //    //}
    //    for (int j=0; j < 3; ++j){
    //         double traction = -this->kb*(2*slmc + 4*mean_curv*(mean_curv*mean_curv 
    //            - gauss_curv))*this->vert_normal[3*ivert+j];
    //        printf("%f  ", traction);
    //    }
    //    printf("\n");
    //}
}

/******************************************************************************/

void asm_set_bc(Assembler* this)
{
    /* Looping over all boundary vertices and setting corresponding */
    /* traction to zero.                                            */
    int num_bverts = bmesh_get_num_vertices(this->bmsh); 
    int ivert;
    for (int ibvert=0; ibvert < num_bverts; ++ibvert){
        ivert = bmesh_get_mapped_entity(this->bmsh, 0, ibvert);
        for (int i=0; i < 3; ++i){
            this->traction[3*ivert+i] = 0.0;
        }
    }
}

/******************************************************************************/

/* Calculates the squared length of all edges */
void asm_calc_edge_len_sq(Assembler* this)
{
    int num_edges = mesh_get_num_entities(this->msh, 1);
    double coords[6];
    double r[3];

    for (int iedge=0; iedge < num_edges; ++iedge){
        mesh_get_ent_coords(this->msh, 1, iedge, coords, NULL);
        for (int i=0; i < 3; ++i){
            r[i] = coords[3*1+i] - coords[3*0+i];
        }
        this->edge_len_sq[iedge] = dot(3, r, r);
    }
}

/******************************************************************************/

/* Calculates the angles of each triangular cell, their cotangents, and  */
/* determines whether a triangle is obtuse. If triangle icell is obtuse, */
/* this->is_obtuse[icell] will contain the local index of the vertex at  */
/* which the angle is obtuse in the cell icell. If triangle icell is not */
/* obtuse this->is_obtuse[icell] is set to -1.                           */
void asm_calc_angles(Assembler* this)
{
    int num_cells = mesh_get_num_entities(this->msh, 2);
    int* conns;
    int nc;
    int is_obtuse;
    double lsq[3];
    double cell_angle[3];
    double cos_cell_angle[3];

    for (int icell=0; icell < num_cells; ++icell){
        /* For triangle cell, nc will return 3 */
        conns = mesh_get_conns(this->msh, 2, 1, icell, &nc);
        /* Using cosine formula */
        for (int j=0; j < nc; ++j){
            lsq[j] = this->edge_len_sq[conns[j]];
        }
        cos_cell_angle[0] = (lsq[1] + lsq[2] - lsq[0])/(2.0*sqrt(lsq[1]*lsq[2]));
        cos_cell_angle[1] = (lsq[2] + lsq[0] - lsq[1])/(2.0*sqrt(lsq[2]*lsq[0]));
        cell_angle[0] = acos(cos_cell_angle[0]);
        cell_angle[1] = acos(cos_cell_angle[1]);
        cell_angle[2] = M_PI - (cell_angle[0] + cell_angle[1]);

        //if (icell == 2) printf("%f %f %f\n", lsq[0], lsq[1], lsq[2]);
        is_obtuse = -1;
        for (int i=0; i < 3; ++i){
            if (cell_angle[i] > M_PI_2){
                is_obtuse = i;
                break;
            }
        }
        this->is_obtuse[icell] = is_obtuse;

        for (int j=0; j < nc; ++j){
            this->cell_angle[3*icell+j] = cell_angle[j];
            this->cell_angle_cot[3*icell+j] = 1.0/tan(cell_angle[j]);
        }
    }
}

/******************************************************************************/

void asm_calc_edge_centers(Assembler* this)
{
    double coords[6];
    int ncoords;
    int nverts;
    int num_edges = mesh_get_num_entities(this->msh, 1);

    for (int iedge=0; iedge < num_edges; ++iedge){
        mesh_get_ent_coords(this->msh, 1, iedge, coords, &ncoords);
        nverts = ncoords/3;
        mat_mean(nverts, 3, coords, 0, this->edge_center+3*iedge);
    }
}

/******************************************************************************/

/* Calculating center of cells -- circumcenter for non-obtuse */
/* triangles, else the midpoint of the side opposite to the   */
/* obtuse angle                                               */
void asm_calc_cell_centers(Assembler* this)
{
    int num_cells = mesh_get_num_entities(this->msh, 2);
    int is_obtuse;
    int g_e;
    double coords[9];

    for (int icell=0; icell < num_cells; ++icell){
        is_obtuse = this->is_obtuse[icell];
    }
    for (int icell=0; icell < num_cells; ++icell){
        is_obtuse = this->is_obtuse[icell];
        if (is_obtuse == -1){
            mesh_get_ent_coords(this->msh, 2, icell, coords, NULL);
            tri_circumcenter(coords, this->cell_center+3*icell,
                    NULL, NULL, NULL);
        }
        else {
            g_e = mesh_get_global_index(this->msh, 2, icell, 1, is_obtuse);
            memcpy(this->cell_center+3*icell, this->edge_center+3*g_e,
                    3*sizeof(double));
        }
    }
}

/******************************************************************************/

/* Calculates the surface gradient of mean curvature */
void asm_calc_sgmc(Assembler* this)
{
    int num_cells = mesh_get_num_entities(this->msh, 2);
    int* conns;
    double coords[9];
    double gxi[3];
    double geta[3];
    double G11, G12, G21, G22;
    double invG11, invG12, invG21, invG22;
    double Hxi, Heta;
    double D;

    for (int icell=0; icell < num_cells; ++icell){
        conns = mesh_get_conns(this->msh, 2, 0, icell, NULL);
        mesh_get_ent_coords(this->msh, 2, icell, coords, NULL);
        /* Calculation of covariant vectors */
        for (int i=0; i < 3; ++i){
            gxi[i] = coords[3*0+i] - coords[3*2+i];
            geta[i] = coords[3*1+i] - coords[3*2+i];
        }

        /* Calculation of metric tensor */
		G11= dot(3, gxi, gxi);
		G12= dot(3, gxi, geta);
		G21= dot(3, geta, gxi);
		G22= dot(3, geta, geta);

		D = G11*G22 - G12*G21;

        /* Calculation of inverse metric tensor */
		invG11 =  G22/D;
		invG12 = -G12/D;
		invG21 = -G21/D;
		invG22 =  G11/D;

        /* surface derivative */
        Hxi = 2.0*this->mean_curv[conns[0]] - 2.0*this->mean_curv[conns[2]];
		Heta = 2.0*this->mean_curv[conns[1]] - 2.0*this->mean_curv[conns[2]]; 
		 
        for (int i=0; i < 3; ++i){
            this->sgmc[3*icell+i] = invG11*Hxi*gxi[i] + invG12*Hxi*geta[i] 
                            + invG21*Heta*gxi[i] + invG22*Heta*geta[i];
         }
    }
}

/******************************************************************************/

// computes the interacial force for a drop */
//void compute_drop_interfacial_force(double fb[beads][3], double curv[beads],
//        double nrm[beads][3])
//{
//  double gamma=KRBC;
//
//  for(int ibead=0; ibead<beads; ibead++){
//    for(int i=0; i<3; ++i){
//      fb[ibead][i] = 2*gamma*curv[ibead]*nrm[ibead][i];
//    }
//  }
//}


/******************************************************************************/
