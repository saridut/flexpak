#include "rotation.h"
#include "triangle_cell.h"
#include "utils_math.h"
#include <string.h>
#include <stdbool.h>

/*****************************************************************************/

double tri_area (const double* coords){

    double a[3];
    double b[3];
    double c[3];

    for (int i=0; i<3; ++i){
        a[i] = coords[3*1+i] - coords[3*0+i];
        b[i] = coords[3*2+i] - coords[3*0+i];
    }
    cross(a, b, c);

    return 0.5*norm(3, c);
}

/*****************************************************************************/

void tri_normal (const double* coords, double* normal){

    double ba[3];
    double ca[3];

    for (int i=0; i<3; ++i){
        ba[i] = coords[3*1+i] - coords[3*0+i];
        ca[i] = coords[3*2+i] - coords[3*0+i];
    }
    cross(ba, ca, normal);

    return unitize(3, normal);
}

/*****************************************************************************/

//From Jonathan Shewchuk's post at
//http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void tri_circumcenter(double* coords,
                        double* cc, double* cr, double* xi, double* eta)
{
    double xba, yba, zba, xca, yca, zca;
    double balength, calength;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;
    double *a, *b, *c;

    a = coords; b = coords+3; c = coords+6; /* Pointers to three vertices */

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;
    
    /* Cross product of these edges. */
    /* Take your chances with floating-point roundoff. */
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;

    /* Calculate the denominator of the formulae. */
    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
                         zcrossbc * zcrossbc);

    /* Calculate offset (from `a') of circumcenter. */
    xcirca = ((balength * yca - calength * yba) * zcrossbc -
              (balength * zca - calength * zba) * ycrossbc) * denominator;
    ycirca = ((balength * zca - calength * zba) * xcrossbc -
              (balength * xca - calength * xba) * zcrossbc) * denominator;
    zcirca = ((balength * xca - calength * xba) * ycrossbc -
              (balength * yca - calength * yba) * xcrossbc) * denominator;

    cc[0] = xcirca;
    cc[1] = ycirca;
    cc[2] = zcirca;

    if (cr != NULL){
        *cr = norm(3, cc);
    }

    /* Get the absolute coordinates of the circumcenter */
    for (int i=0; i < 3; ++i){
        cc[i] += a[i];
    }

    if (xi != (double *) NULL) {
      /* To interpolate a linear function at the circumcenter, define a     */
      /*   coordinate system with a xi-axis directed from `a' to `b' and    */
      /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
      /*   are computed by Cramer's Rule for solving systems of linear      */
      /*   equations.                                                       */

      /* There are three ways to do this calculation - using xcrossbc, */
      /*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
      /*   magnitude, to improve stability and avoid division by zero. */
      if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
          ((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc))) {
        *xi = (ycirca * zca - zcirca * yca) / xcrossbc;
        *eta = (zcirca * yba - ycirca * zba) / xcrossbc;
      } else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc)) {
        *xi = (zcirca * xca - xcirca * zca) / ycrossbc;
        *eta = (xcirca * zba - zcirca * xba) / ycrossbc;
      } else {
        *xi = (xcirca * yca - ycirca * xca) / zcrossbc;
        *eta = (ycirca * xba - xcirca * yba) / zcrossbc;
      }
    }
}

/*****************************************************************************/

/* Calculates the projection operator `op' for a given normal. */
/* op = I - pp, where p is the normal vector.                  */
void op_proj(const double* normal, double* op){

    outer(3, 3, normal, normal, op);

    /*setting op = -op */
    for (int i=0; i<3; ++i){
        for (int j=0; i<3; ++j){
            op[3*i+j] = -op[3*i+j];
        }
    }

    /* Adding one to the diagonal elements, i.e. op[i,i] = 1 - op[i,i] */
    for (int i=0; i<3; ++i){
        op[3*i+i] += 1.0;
    }
}

/*****************************************************************************/

/* Calculates coordinates of vertices in local basis */
void tri_to_local(const double* g_verts, double* l_verts, double* dcm)
{
    double tmp_verts[9];
    double l_basis[9];
    double dx[3];

    /* Translate such that vertex 0 is at the global origin */
    memcpy(tmp_verts, g_verts, 9*sizeof(double));
    memcpy(dx, &g_verts[0], 3*sizeof(double));
    for (int i=0; i < 3; ++i){
        dx[i] = -dx[i];
    }
    translate(3, tmp_verts, &dx[0]);

    tri_get_local_basis(tmp_verts, l_basis);
    /* For two orthogonal basis sets A and B, one should next call  */
    /* dcm_from_axes(A, B, dcm) to obtain the dcm of B w.r.t. A. In */
    /* case A is identity (which is the case here), dcm = B, so no  */
    /* function call is necessary.                                  */
    shift_vectors_dcm(3, tmp_verts, l_basis, false, l_verts);

    if (dcm != NULL){
        memcpy(dcm, l_basis, 9*sizeof(double));
    }
}

/*****************************************************************************/

void tri_get_local_basis(const double* verts, double* basis)
{
    double v0[3];
    double v1[3];
    double v2[3];

    for (int i=0; i < 3; ++i){
        v0[i] = verts[0+i];
        v1[i] = verts[3+i];
        v2[i] = verts[6+i];
    }

    double* xhat = &basis[0];
    double* yhat = &basis[3];
    double* zhat = &basis[6];

    for (int j=0; j < 3; ++j){
        xhat[j] = v1[j]- v0[j];
        yhat[j] = v2[j]- v1[j];
    }
    cross(xhat, yhat, zhat);
    cross(zhat, xhat, yhat);
    unitize(3, xhat);
    unitize(3, yhat);
    unitize(3, zhat);

}
/*****************************************************************************/

/* Calculates the displacement of triangle B w.r.t. triangle A, both */
/* expressed in the same basis.                                      */
void tri_get_displacement(const double* verts_A, const double* verts_B,
        double* u, double* v)
{
    for (int i=0; i < 3; ++i){
        u[i] = verts_B[3*i] - verts_A[3*i];
        v[i] = verts_B[3*i+1] - verts_A[3*i+1];
    }
}


/*****************************************************************************/
