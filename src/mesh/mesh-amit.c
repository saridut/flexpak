#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <time.h>
#include <fftw3.h>
#include "globalvars.h"
#include "sharedvars.h"
#include "prototype.h"
#include "mesh.h"

/*****************************************************************************/

void init_normal_beads (double normal[NTRIANGLE_MAX][3], double nrm[beads][3], int ipart)
{
    double wx[beads], wy[beads], wz[beads];
    int i;
    int i1,i2,i3;
    int itriangle;
    double norm;
    
    /* Initialize weights for the global solution */
    for(i=0;i<beads;i++) {
        wx[i]=0.0;
        wy[i]=0.0;
        wz[i]=0.0;
    }
    
    /* Store weights for the global calculation */
    /* Also used for computing the area weighted normal average for each bead */
    for(itriangle=0;itriangle<ntriangles;itriangle++) {
        i1 = triangle[itriangle][0];
        i2 = triangle[itriangle][1];
        i3 = triangle[itriangle][2];
        
        wx[i1] += area[ipart][itriangle]/3.0*normal[itriangle][0];
        wx[i2] += area[ipart][itriangle]/3.0*normal[itriangle][0];
        wx[i3] += area[ipart][itriangle]/3.0*normal[itriangle][0];
        
        wy[i1] += area[ipart][itriangle]/3.0*normal[itriangle][1];
        wy[i2] += area[ipart][itriangle]/3.0*normal[itriangle][1];
        wy[i3] += area[ipart][itriangle]/3.0*normal[itriangle][1];
        
        wz[i1] += area[ipart][itriangle]/3.0*normal[itriangle][2];
        wz[i2] += area[ipart][itriangle]/3.0*normal[itriangle][2];
        wz[i3] += area[ipart][itriangle]/3.0*normal[itriangle][2];
    
        normal_t[ipart][itriangle][0] = normal[itriangle][0];
        normal_t[ipart][itriangle][1] = normal[itriangle][1];
        normal_t[ipart][itriangle][2] = normal[itriangle][2];
    }
    
    /*Computing the area weighted normal average for each bead */
    for(i=0;i<beads;i++) {
        norm = wx[i]*wx[i] + wy[i]*wy[i] + wz[i]*wz[i];
        norm = sqrt(norm);
        
        nrm[i][0] = wx[i]/norm;
        nrm[i][1] = wy[i]/norm;
        nrm[i][2] = wz[i]/norm;
    }
}

/*****************************************************************************/

void compute_normal_curvature (double nrm[beads][3], double curv[beads],
        double xb[beads][3], int beads_to_beads[beads][6])
{
  int i,j,k;
  int ibead,jbead;
  int i1,i2,i3;
  int itriangle;
  double norm;
  double xneigh[6][3];
  double x_tr[6][3];
  double Ct[3][3];
  double u[3],Q[3][3];
  double xp[3],yp[3],zp[3];
  double A1[5][5],brhs1[5],rdist1[5];
  double A2[6][6],brhs2[6],rdist2[6];
  int INFO, N, NRHS, IPIV[6],LDA, LDB, M, RANK;
  int LWORK=50;
  double WORK[50];
  double RCOND=1E-10;
  double S[5];
  double nrm_l[3],nrm_g[3];
  double error;
  FILE *fp;

  /*fp = fopen("nrm.txt","w+");
  for(ibead=0;ibead<beads;ibead++)
  {
    for(i=0;i<3;i++)
    {
      fprintf(fp,"%E\t",nrm[ibead][i]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);*/

  /* First for the first 12 vertices with coordination number 5 */
  for(ibead=0;ibead<12;ibead++)
  {
    /* store coordinates of the connected beads relative to ibead */
    for(i=0;i<5;i++)
    {
      jbead = beads_to_beads[ibead][i];
      for(j=0;j<3;j++)
      {
        xneigh[i][j] = xb[jbead][j]-xb[ibead][j];
      }
    }

    do
    {
      /* Store local z axis coordinates in global frame */
      zp[0] = nrm[ibead][0];
      zp[1] = nrm[ibead][1];
      zp[2] = nrm[ibead][2];

      /* find the rotation matrix: Use householder's algorithm */
      if (fabs(zp[2]-1) > 1E-10)
      {
        u[0] = 0.0 - zp[0];
        u[1] = 0.0 - zp[1];
        u[2] = 1.0 - zp[2];

        norm = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        
        /* normalize */
        u[0] = u[0]/norm;
        u[1] = u[1]/norm;
        u[2] = u[2]/norm;

        /* Find the Q matrix */
        for(i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            if (i==j)
            {
              Q[i][j] = 1.0 - 2.0*u[i]*u[j];
            }
            else
            {
              Q[i][j] = -2.0*u[i]*u[j];
            }
          }
        }
      }
      else // no rotation necessary Q=I
      {
        
        for(i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            if (i==j)
            {
              Q[i][j] = 1.0;
            }
            else
            {
              Q[i][j] = 0.0;
            }
          }
        }
      }


      /* find coordinates of the local x axis in global system = Q^T*(1,0,0) */
      xp[0] = Q[0][0];
      xp[1] = Q[1][0];
      xp[2] = Q[2][0];

      /* find yp = z * x*/
      cross_product(zp,xp,yp);

      /* compute the coordinate transformation matrix: C^T */
      for(i=0;i<3;i++)
      {
        Ct[0][i] = xp[i];
        Ct[1][i] = yp[i];
        Ct[2][i] = zp[i];
      }

      /* transform coordinates of neighboring beads to local coordinates */
      for(i=0;i<5;i++)
      {
        for(j=0;j<3;j++)
        {
          x_tr[i][j] = 0.0;
          for(k=0;k<3;k++)
          {
            x_tr[i][j] += Ct[j][k]*xneigh[i][k];
          }
        }
      }

      /* Find the distance from the center */
      for(i=0;i<5;i++)
      {
        rdist1[i] = sqrt(x_tr[i][0]*x_tr[i][0] + x_tr[i][1]*x_tr[i][1] + x_tr[i][2]*x_tr[i][2]);
      }

      /* build the A matrix, store in fortran column major format */
      for(i=0;i<5;i++)
      {
        A1[0][i] = x_tr[i][0]/rdist1[i];
        A1[1][i] = x_tr[i][1]/rdist1[i];
        A1[2][i] = x_tr[i][0]*x_tr[i][0]/rdist1[i];
        A1[3][i] = x_tr[i][0]*x_tr[i][1]/rdist1[i];
        A1[4][i] = x_tr[i][1]*x_tr[i][1]/rdist1[i];

        brhs1[i] = x_tr[i][2]/rdist1[i];
      }

      /* Solve using dgesv */
      N=5;NRHS=1;LDA=5;LDB=5;
      dgesv_(&N,&NRHS,A1,&LDA,IPIV,brhs1,&LDB,&INFO);

      /* find new normal in local coordinates */
      norm = sqrt(1.0 + brhs1[0]*brhs1[0] + brhs1[1]*brhs1[1]);
      nrm_l[0] = -brhs1[0]/norm;
      nrm_l[1] = -brhs1[1]/norm;
      nrm_l[2] = 1.0/norm;

      /* transform to global coordinates */
      for(i=0;i<3;i++)
      {
        nrm_g[i]=0.0;
        for(j=0;j<3;j++)
        {
          nrm_g[i] += Ct[j][i]*nrm_l[j];
        }
      }

      /* find error */
      error=0.0;
      for(i=0;i<3;i++)
      {
        error += pow(zp[i]-nrm_g[i],2);
      }
      error = sqrt(error);

      /* store in normal array */
      for(i=0;i<3;i++)
      {
        nrm[ibead][i] = nrm_g[i];
      }


    } while (error > 1E-3);

    /* store curvature */
    curv[ibead] = -brhs1[2] - brhs1[4];


//     printf("%d\t%E\n",ibead,curv[ibead]);


  }

  /* Remaining vertices with coordination number 6 */
  for(ibead=12;ibead<beads;ibead++)
  {
    /* store coordinates of the connected beads relative to ibead */
    for(i=0;i<6;i++)
    {
      jbead = beads_to_beads[ibead][i];
      for(j=0;j<3;j++)
      {
        xneigh[i][j] = xb[jbead][j]-xb[ibead][j];
      }
    }

    do
    {
      /* Store zp axis coordinates */
      zp[0] = nrm[ibead][0];
      zp[1] = nrm[ibead][1];
      zp[2] = nrm[ibead][2];

      
      if (fabs(zp[2]-1) > 1E-10)
      {
        u[0] = 0.0 - zp[0];
        u[1] = 0.0 - zp[1];
        u[2] = 1.0 - zp[2];
        
        norm = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        
        /* normalize */
        u[0] = u[0]/norm;
        u[1] = u[1]/norm;
        u[2] = u[2]/norm;
        
        /* Find the Q matrix */
        for(i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            if (i==j)
            {
              Q[i][j] = 1.0 - 2.0*u[i]*u[j];
            }
            else
            {
              Q[i][j] = -2.0*u[i]*u[j];
            }
          }
        }
      }
      else // no rotation necessary Q=I
      {
        
        for(i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            if (i==j)
            {
              Q[i][j] = 1.0;
            }
            else
            {
              Q[i][j] = 0.0;
            }
          }
        }
      }

      /* find coordinates of xp axis in global system */
      xp[0] = Q[0][0];
      xp[1] = Q[1][0];
      xp[2] = Q[2][0];

      /* find yp = z * x*/
      cross_product(zp,xp,yp);

      /* compute the coordinate transformation matrix: C^T */
      for(i=0;i<3;i++)
      {
        Ct[0][i] = xp[i];
        Ct[1][i] = yp[i];
        Ct[2][i] = zp[i];
      }

      /* transform coordinates of neighboring beads to local coordinates */
      for(i=0;i<6;i++)
      {
        for(j=0;j<3;j++)
        {
          x_tr[i][j] = 0.0;
          for(k=0;k<3;k++)
          {
            x_tr[i][j] += Ct[j][k]*xneigh[i][k];
          }
        }
      }

      /* Find the distance from the center */
      for(i=0;i<6;i++)
      {
        rdist2[i] = sqrt(x_tr[i][0]*x_tr[i][0] + x_tr[i][1]*x_tr[i][1] + x_tr[i][2]*x_tr[i][2]);
      }

      /* build the A matrix, store in fortran column major format */
      for(i=0;i<6;i++)
      {
        A2[0][i] = x_tr[i][0]/rdist2[i];
        A2[1][i] = x_tr[i][1]/rdist2[i];
        A2[2][i] = x_tr[i][0]*x_tr[i][0]/rdist2[i];
        A2[3][i] = x_tr[i][0]*x_tr[i][1]/rdist2[i];
        A2[4][i] = x_tr[i][1]*x_tr[i][1]/rdist2[i];

        brhs2[i] = x_tr[i][2]/rdist2[i];
      }

      /* Solve using dgelss */
      M=6;N=5;NRHS=1;LDA=M;LDB=M;
      dgelss_(&M,&N,&NRHS,A2,&LDA,brhs2,&LDB,S,&RCOND,&RANK,WORK,&LWORK,&INFO);


      /* find new normal in local coordinates */
      norm = sqrt(1.0 + brhs2[0]*brhs2[0] + brhs2[1]*brhs2[1]);
      nrm_l[0] = -brhs2[0]/norm;
      nrm_l[1] = -brhs2[1]/norm;
      nrm_l[2] = 1.0/norm;

      /* transform to global coordinates */
      for(i=0;i<3;i++)
      {
        nrm_g[i]=0.0;
        for(j=0;j<3;j++)
        {
          nrm_g[i] += Ct[j][i]*nrm_l[j];
        }
      }

      /* find error */
      error=0.0;
      for(i=0;i<3;i++)
      {
        error += pow(zp[i]-nrm_g[i],2);
      }
      error = sqrt(error);

      /* store in normal array */
      for(i=0;i<3;i++)
      {
        nrm[ibead][i] = nrm_g[i];
      }


    } while (error > 1E-3);


    /* store curvature */
    curv[ibead] = -brhs2[2] - brhs2[4];
//     printf("%d\t%E\n",ibead,curv[ibead]);


  }

/*fp = fopen("nrm1.txt","w+");
for(ibead=0;ibead<beads;ibead++)
{
  for(i=0;i<3;i++)
  {
    fprintf(fp,"%E\t",nrm[ibead][i]);
  }
  fprintf(fp,"\n");
}
fclose(fp);
exit(0);*/

}

/*****************************************************************************/

/*-------------------------------------------------------------
 * Computes the solid angle subtended by a triangle at
 * one of its own vertex  (see ieee transactions of biomedical
 * engineering, 45, 980, 1998)
 - *------------------------------------------------------------*/
void solid_angle(double *Omega, double x1[3],double x2[3], double x3[3], double n[3])
{
    int i,j,k;
    double pi;
    double x[3],y[3];
    double z[3];
    double ndotx,ndoty,ndotz,xdoty;
    double xm,ym;
    double nr,dr;
    
    pi = 4.0*atan(1.0);
    
    /* find sides of the triangle relative to the origin vertex */
    for(i=0;i<3;i++) {
        x[i]=x2[i]-x1[i];
        y[i]=x3[i]-x1[i];
    }
    
    cross_product(x,y,z);
    
    xm  = sqrt(dot_product(3,x,x));
    ym  = sqrt(dot_product(3,y,y));
    
    ndotx = dot_product(3,n,x);
    ndoty = dot_product(3,n,y);
    ndotz = dot_product(3,n,z);
    xdoty = dot_product(3,x,y);

    nr  = -2.0*ndotz*(ndotx*ym + ndoty*xm);
    dr = pow(xm*ym + xdoty,2) - pow(ndotx*ym + ndoty*xm,2) + pow(ndotz,2);
    *Omega = atan(nr/dr);
}

/******************************************************************************/

/*-----------------------------------------------------------------------------------
uin = coefficient of the Green's function
uout = output of the double layer integral
-------------------------------------------------------------------------------------*/
void prvec_double_layer(double xb[npart][beads][3],double uin[npart][beads][3],
                        double nrm[npart][beads][3], double ub[npart][beads][3])
{

  int i,j,k,ibead;
  int itriangle,jtriangle;
  double  xi[3];
  int i1,i2,i3;
  double u[3],u1[3];
  double w[npart][beads];
  double uxf[nx][ny][nz], uyf[nx][ny][nz], uzf[nx][ny][nz],pf[nx][ny][nz];
  double ubc1[nx][ny][3], ubc2[nx][ny][3];
  double gaussx[nx][ny][nz],gaussy[nx][ny][nz],gaussz[nx][ny][nz],gaussz_p[nx][ny][nz];
  FILE *fp;
  char name[20];
  double fb[beads][3];
  double sigma_g[nx][ny][nz][3][3];
  double pi;
  double uint[beads][3];
  int ix,iy,iz;
  double dx,dy;
  int member;
  int ipart,jpart;
  double Omega;
	double dudx[nx][ny][nz],dudy[nx][ny][nz],dudz[nx][ny][nz];
	double dvdx[nx][ny][nz],dvdy[nx][ny][nz],dvdz[nx][ny][nz];
	double dwdx[nx][ny][nz],dwdy[nx][ny][nz],dwdz[nx][ny][nz];
	
	
  
  
  pi = 4.0*atan(1.0);
  
  /* Mesh size */
  dx = Lx/nx;
  dy = Ly/ny;

  /* Initialize velocity to zero */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<beads;i++)
    {
      for(j=0;j<3;j++)
      {
        ub[ipart][i][j]=0.0;
      }
    }
  }
  
  /* Initialize weights for the global solution */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<beads;i++)
    {
      w[ipart][i]=0.0;
    }
  }


  /* Store weights for the global calculation */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];
      
      w[ipart][i1] += area[ipart][itriangle]/3.0;
      w[ipart][i2] += area[ipart][itriangle]/3.0;
      w[ipart][i3] += area[ipart][itriangle]/3.0;
      
    }
  }

#if (OMP == 1)
  omp_set_num_threads(NTHREADS);
#endif


  /* Sum contribution from the self element at each of the nodes*/
  /* Also find the weight for each node for the global solution computaiton */
#if (npart > 1)
#if (OMP == 1)
#pragma omp parallel for private(ipart,jpart,ibead,j,itriangle,jtriangle,xi,i1,i2,i3,member,u,Omega)
#endif
#endif
  for(ipart=0;ipart<npart;ipart++)
  {
#if (npart == 1)
#if (OMP == 1)
#pragma omp parallel for private(jpart,ibead,j,itriangle,jtriangle,xi,i1,i2,i3,member,u,Omega)
#endif
#endif
    for(ibead=0;ibead<beads;ibead++)
    {
      /* Store position of the field point */
      for(i=0;i<3;i++)
      {
        xi[i] = xb[ipart][ibead][i];
      }
      
      /* compute velocity at this point due to all triangular elements */
      for(itriangle=0;itriangle<triangle_count[ipart][ibead];itriangle++)
      {
        jpart = triangle_list[ipart][ibead][itriangle][0];
        jtriangle = triangle_list[ipart][ibead][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0;
        
        //       printf("%lf\n",area[triangle_list[ibead][itriangle]]);
        
        if (ibead == i1 && jpart == ipart)
        {
          member=1;
        }
        else if (ibead == i2  && jpart == ipart)
        {
          i2=i1;
          i1 = ibead;
          member=1;
        }
        else if (ibead == i3  && jpart == ipart)
        {
          i3=i1;
          i1=ibead;
          member=1;
        }
        
        /* Perform the integration */
        if (member==0)
        {
          integrate5(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],uin[jpart][i3],
                     u,alpha,nrm[ipart][ibead]);
        }
        else
        {
          if (SOLID_ANGLE == 1)
          {
            /* compute solid angle */
            solid_angle(&Omega,xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],nrm[ipart][ibead]);
          
            integrate5c(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],
                        uin[jpart][i2],uin[jpart][i3],u,alpha,normal_t[ipart][jtriangle],Omega);
          }
          else
          {
          /* using normal of the triangle as the normal of ibead, this makes it non-singular */
            integrate5a(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],
                      uin[jpart][i3],u,alpha,nrm[ipart][ibead],normal_t[ipart][jtriangle]);
          }
          
        }
        
        
        /* Sum to total velocity */
        for(i=0;i<3;i++)
        {
          ub[ipart][ibead][i] += u[i];
        }
      }
      
    }
  }
  




/*----------------------- Global Solution Calculation ------------------------------*/

  /*----- Boundary Conditions ---- */
  
  /* Initialize */
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<3;k++)
      {
        ubc1[i][j][k]=0.0;
        ubc2[i][j][k]=0.0;
      }
    }
  }


  /* set ubc = - u_l for poiseuille flow */
  /* Top plate */
  #if (OMP == 1)
  #pragma omp parallel for private(ix,iy,xi,itriangle,jpart,jtriangle,i1,i2,i3,member,u)
  #endif
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = Lz;

      for(itriangle=0;itriangle<triangle_count_t[ix][iy];itriangle++)
      {
        jpart = triangle_list_t[ix][iy][itriangle][0];
        jtriangle = triangle_list_t[ix][iy][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];

        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],uin[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],uin[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc2[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }

  /* Bottom plate */
  #if (OMP == 1)
  #pragma omp parallel for private(ix,iy,xi,itriangle,jpart,jtriangle,i1,i2,i3,member,u)
  #endif
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = 0.0;

      for(itriangle=0;itriangle<triangle_count_b[ix][iy];itriangle++)
      {
        jpart = triangle_list_b[ix][iy][itriangle][0];
        jtriangle = triangle_list_b[ix][iy][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];

        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],uin[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,uin[jpart][i1],uin[jpart][i2],uin[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc1[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }


  /* Distribute density to mesh */
  /* Initialize mesh force density */
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      for(iz=0;iz<nz;iz++)
      {
        gaussx[ix][iy][iz] = 0.0;
        gaussy[ix][iy][iz] = 0.0;
        gaussz[ix][iy][iz] = 0.0;
				gaussz_p[ix][iy][iz] = 0.0;
      }
    }
  }
 

  for(ipart=0;ipart<npart;ipart++)
  {
    distribute_density(xb[ipart],uin[ipart],gaussx,gaussy,gaussz,gaussz_p,w[ipart],alpha,xz);
  }
  
 
 
	/* solve for the velocity and pressure at the mesh points */
	global_velocity_inhomogeneous(gaussx,gaussy,gaussz,gaussz_p,ubc1,ubc2,uxf,uyf,uzf,pf,
																dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,															
																xxi,cz,
																u1H,v1H,w1H,dwdz1H,dudz1H,dvdz1H,p1H,
															  u2H,v2H,w2H,dwdz2H,dudz2H,dvdz2H,p2H);
	

  /* Compute the global stress tensor at the mesh points */
	
#if (STRESS_HI == 1)
	compute_global_stress_tensor_spectral(dudx,dudy,dudz,dvdx,dvdy,dvdz,
																	dwdx,dwdy,dwdz,pf,sigma_g);
#else
  compute_global_stress_tensor(uxf,uyf,uzf,pf,sigma_g,xz);
#endif


  /* Store first column of sigma in uxf, uyf, uzf to interpolate */
  for (i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
        uxf[i][j][k] = sigma_g[i][j][k][0][0];
        uyf[i][j][k] = sigma_g[i][j][k][1][0];
        uzf[i][j][k] = sigma_g[i][j][k][2][0];
      }
    }
  }

/* Interpolate velocity from the mesh to the bead position: store in uint */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        uint[ibead][i]=0.0;
      }
    }
    
#if (INTERP == 2)
    interpolate2(xb[ipart],uint,uxf,uyf,uzf,xz);
#else
		interpolate4(xb[ipart],uint,uxf,uyf,uzf,xz);
#endif
		
    /* multiply by nrm[ibead][0] */
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        ub[ipart][ibead][i] += uint[ibead][i]*nrm[ipart][ibead][0];
      }
    }
  }

  /* Store second column of sigma in uxf, uyf, uzf to interpolate */
  for (i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
        uxf[i][j][k] = sigma_g[i][j][k][0][1];
        uyf[i][j][k] = sigma_g[i][j][k][1][1];
        uzf[i][j][k] = sigma_g[i][j][k][2][1];
      }
    }
  }


  /* Interpolate velocity from the mesh to the bead position: store in uint */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        uint[ibead][i]=0.0;
      }
    }
    
#if (INTERP == 2)
    interpolate2(xb[ipart],uint,uxf,uyf,uzf,xz);
#else
		interpolate4(xb[ipart],uint,uxf,uyf,uzf,xz);
#endif
		
    /* multiply by nrm[ibead][1] */
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        ub[ipart][ibead][i] += uint[ibead][i]*nrm[ipart][ibead][1];
      }
    }
  }

  /* Store third column of sigma in uxf, uyf, uzf to interpolate */
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<nz;k++)
      {
        uxf[i][j][k] = sigma_g[i][j][k][0][2];
        uyf[i][j][k] = sigma_g[i][j][k][1][2];
        uzf[i][j][k] = sigma_g[i][j][k][2][2];
      }
    }
  }

  /* Interpolate velocity from the mesh to the bead position: store in uint */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        uint[ibead][i]=0.0;
      }
    }
#if (INTERP == 2)
    interpolate2(xb[ipart],uint,uxf,uyf,uzf,xz);
#else
		interpolate4(xb[ipart],uint,uxf,uyf,uzf,xz);
#endif
    /* multiply by nrm[ibead][2] */
    for(ibead=0;ibead<beads;ibead++)
    {
      for(i=0;i<3;i++)
      {
        ub[ipart][ibead][i] += uint[ibead][i]*nrm[ipart][ibead][2];
      }
    }
  }

}

/******************************************************************************/

void prvec(double xb[npart][beads][3], double fb[npart][beads][3], double ub[npart][beads][3],
	double nrm_b[npart][beads][3])
{
  
    int i,j,k;
    int ibead,itriangle,jtriangle;
    double  xi[3];
    double r1;
    int i1,i2,i3;
    double pi;
    double u[3];
    int member;
    double weights[npart][beads];
    double uxf[nx][ny][nz], uyf[nx][ny][nz], uzf[nx][ny][nz], pf[nx][ny][nz];
    double ubc1[nx][ny][3], ubc2[nx][ny][3];
    double gaussx[nx][ny][nz],gaussy[nx][ny][nz],gaussz[nx][ny][nz],gaussz_p[nx][ny][nz];
	double dudx[nx][ny][nz],dudy[nx][ny][nz],dudz[nx][ny][nz];
	double dvdx[nx][ny][nz],dvdy[nx][ny][nz],dvdz[nx][ny][nz];
	double dwdx[nx][ny][nz],dwdy[nx][ny][nz],dwdz[nx][ny][nz];
    FILE *fp;
    char name[20];
    int ix,iy,iz;
    double dx,dy;
    int ipart,jpart;
	double fn,f1[3],f2[3],f3[3];
	
  pi = 4.0*atan(1.0);

  /* Mesh size */
  dx = Lx/nx;
  dy = Ly/ny;
  
  
  /* Initialize velocity to zero */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<beads;i++)
    {
      for(j=0;j<3;j++)
      {
        ub[ipart][i][j]=0.0;
      }
    }
  }

#if (OMP == 1)
omp_set_num_threads(NTHREADS);
#endif


  
  /* Contribution from local force density */
#if (npart > 1)
#if (OMP == 1)
#pragma omp parallel for private(ipart,jpart,ibead,j,itriangle,jtriangle,xi,i1,i2,i3,member,u)
#endif
#endif
  for(ipart=0;ipart<npart;ipart++)
  {
#if (npart == 1)
#if (OMP == 1)
#pragma omp parallel for private(jpart,ibead,j,itriangle,jtriangle,xi,i1,i2,i3,member,u)
#endif
#endif
    for(ibead=0;ibead<beads;ibead++)
    {
      for(j=0;j<3;j++)
      {
        xi[j] = xb[ipart][ibead][j];
      }
      
#if (SING_SUBT == 1)
      fn = fb[ipart][ibead][0]*nrm_b[ipart][ibead][0] +
					 fb[ipart][ibead][1]*nrm_b[ipart][ibead][1] +
					 fb[ipart][ibead][2]*nrm_b[ipart][ibead][2];
#endif
      
      for(itriangle=0;itriangle<triangle_count[ipart][ibead];itriangle++)
      {
        jpart = triangle_list[ipart][ibead][itriangle][0];
        jtriangle = triangle_list[ipart][ibead][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0;
        
        //       printf("%lf\n",area[triangle_list[ibead][itriangle]]);
        
        if (ibead == i1 && ipart == jpart)
        {
          member=1;
        }
        else if (ibead == i2 && jpart == ipart)
        {
          i2=i1;
          i1 = ibead;
          member=1;
        }
        else if (ibead == i3 && jpart == ipart)
        {
          i3=i1;
          i1=ibead;
          member=1;
        }
        
        
#if (SING_SUBT == 0 )
        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
#else
				for(j=0;j<3;j++)
				{
					f1[j] = fb[jpart][i1][j] - fn*normal_t[jpart][jtriangle][j];
					f2[j] = fb[jpart][i2][j] - fn*normal_t[jpart][jtriangle][j];
					f3[j] = fb[jpart][i3][j] - fn*normal_t[jpart][jtriangle][j];
				}
				
        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,f1,f2,f3,u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,f1,f2,f3,u,alpha);
        }
#endif
        
        for(i=0;i<3;i++)
        {
          ub[ipart][ibead][i] += u[i];
        }
        
      }
    }
  }
  

/*	fp = fopen("U.txt","w+");
	for(ipart=0;ipart<npart;ipart++)
	{
		for(ibead=0;ibead<beads;ibead++)
		{
			fprintf(fp,"%d\t%E\t%E\t%E\n",ibead,ub[ipart][ibead][0],ub[ipart][ibead][0],ub[ipart][ibead][0]);
		}
	}
	fclose(fp);
	exit(0);*/

//      return;
  
  /*---------------- Global solution calculation -----------------*/


  /*----- Boundary Conditions ---- */
  /* Initialize */
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<3;k++)
      {
        ubc1[i][j][k]=0.0;
        ubc2[i][j][k]=0.0;
      }
    }
  }


  /* set ubc = - u_l for poiseuille flow */
  /* Top plate */
  #if (OMP == 1)
  #pragma omp parallel for private(ix,iy,xi,itriangle,jpart,jtriangle,i1,i2,i3,member,u)
  #endif
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = Lz;

      for(itriangle=0;itriangle<triangle_count_t[ix][iy];itriangle++)
      {
        jpart = triangle_list_t[ix][iy][itriangle][0];
        jtriangle = triangle_list_t[ix][iy][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc2[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }

  /* bottom plate */
#if (OMP == 1)
#pragma omp parallel for private(ix,iy,xi,itriangle,jpart,jtriangle,i1,i2,i3,member,u)
#endif
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = 0.0;

      for(itriangle=0;itriangle<triangle_count_b[ix][iy];itriangle++)
      {
        jpart = triangle_list_b[ix][iy][itriangle][0];
        jtriangle = triangle_list_b[ix][iy][itriangle][1];
        
        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc1[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }
  
  /* Find weights for the global solution */
  /* Initialize weights for the global solution */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<beads;i++)
    {
      weights[ipart][i]=0.0;
    }
  }
  for(ipart=0;ipart<npart;ipart++)
  {
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];
      
      /* Store weights for the global calculation */
      weights[ipart][i1] += area[ipart][itriangle]/3.0;
      weights[ipart][i2] += area[ipart][itriangle]/3.0;
      weights[ipart][i3] += area[ipart][itriangle]/3.0;
      
    }
  }

  /* Distribute density to mesh */
  /* Initialize mesh force density */
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      for(iz=0;iz<nz;iz++)
      {
        gaussx[ix][iy][iz] = 0.0;
        gaussy[ix][iy][iz] = 0.0;
        gaussz[ix][iy][iz] = 0.0;
				gaussz_p[ix][iy][iz] = 0.0;
      }
    }
  }

  for(ipart=0;ipart<npart;ipart++)
  {
    distribute_density(xb[ipart],fb[ipart],gaussx,gaussy,gaussz,gaussz_p,weights[ipart],alpha,xz);
  }
  
 /* printf("Here\n");
  fp =fopen("rho.txt","w+");
	for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      for(iz=0;iz<nz;iz++)
      {
        fprintf(fp,"%E\t%E\t%E\t%E\n",gaussx[ix][iy][iz],gaussy[ix][iy][iz],gaussz[ix][iy][iz],gaussz_p[ix][iy][iz]);
      }
    }
  }
	fclose(fp); */
  
  /* solve for the velocity and pressure at the mesh points */
	global_velocity_inhomogeneous(gaussx,gaussy,gaussz,gaussz_p,ubc1,ubc2,uxf,uyf,uzf,pf,
																dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,															
																xxi,cz,
																u1H,v1H,w1H,dwdz1H,dudz1H,dvdz1H,p1H,
															  u2H,v2H,w2H,dwdz2H,dudz2H,dvdz2H,p2H);

  
	
 /* fp =fopen("Vel.txt","w+");
	for(iz=0;iz<nz;iz++)
  {
		fprintf(fp,"%E\t%E\t%E\t%E\n",xz[iz],uxf[nx/3][ny/3][iz],uyf[nx/3][ny/3][iz],uzf[nx/3][ny/3][iz]);
  }
  fclose(fp);
  exit(0); */
	
  /* Interpolate velocity from the mesh to the bead position */
  for(ipart=0;ipart<npart;ipart++)
  {
#if (INTERP == 2)
    interpolate2(xb[ipart],ub[ipart],uxf,uyf,uzf,xz);
#else
		interpolate4(xb[ipart],ub[ipart],uxf,uyf,uzf,xz);
#endif
  }
  
  /*fp = fopen("U.txt","w+");
	for(ipart=0;ipart<npart;ipart++)
	{
		for(ibead=0;ibead<beads;ibead++)
		{
			fprintf(fp,"%d\t%E\t%E\t%E\n",ibead,ub[ipart][ibead][0],ub[ipart][ibead][1],ub[ipart][ibead][2]);
		}
	}
	fclose(fp);
	exit(0);*/
  
  
}

/*------------------------------------------------------------------
 Find the particle pairs which overlap or violate the minimum
 gap specification
 -------------------------------------------------------------------*/ 
void find_overlaps_gap(double xb[npart][beads][3], double nrm_b[npart][beads][3],int overlap_pairs[npart][npart],
	double min_gap[npart][npart],int *overlaps)
{
	int ipart,jpart;
	int ibead,jbead;
	int i,j;
	double xi[3],xj[3];
	double x,y,z,r;
	double xc,yc,zc;
	double dotp;	

	*overlaps=0;
	
	for(ipart=0;ipart<npart;ipart++)
	{
		for(jpart=0;jpart<npart;jpart++)
		{
			overlap_pairs[ipart][jpart]=0;
			min_gap[ipart][jpart]=GAP_MIN;
		}
	}

	for(i=0;i<ovp_pair_count;i++)
	{
		ipart = ovp_list[i][0];
		ibead = ovp_list[i][1];
		jpart = ovp_list[i][2];
		jbead = ovp_list[i][3];
		
		for(j=0;j<3;j++)
		{
			xi[j] = xb[ipart][ibead][j];
			xj[j] = xb[jpart][jbead][j];
		}
		
		// difference vector
		x = xi[0]-xj[0];
		y = xi[1]-xj[1];
		z = xi[2]-xj[2];
		
		/* correct for periodicity */
		x  = x - Lx*floor(x/Lx+0.5);
        y  = y - Ly*floor(y/Ly+0.5);
		r = sqrt(x*x+y*y+z*z); 
		
		/* find dot_product to check overlap */
		dotp = (x*nrm_b[jpart][jbead][0] + y*nrm_b[jpart][jbead][1] + z*nrm_b[jpart][jbead][2])/r;
		
		if (dotp < 0.0) // overlaps
		{
			overlap_pairs[ipart][jpart]=1;
			overlap_pairs[jpart][ipart]=1;
			
			if ( -r < min_gap[ipart][jpart])
			{
				min_gap[ipart][jpart]=-r;
				min_gap[jpart][ipart]=-r;			
			}
			
			*overlaps = *overlaps + 1;
		}			
		else if (r < GAP_MIN)
		{
			overlap_pairs[ipart][jpart]=1;
			overlap_pairs[jpart][ipart]=1;
			
			if ( r < min_gap[ipart][jpart])
			{
				min_gap[ipart][jpart]=r;
				min_gap[jpart][ipart]=r;			
			}	
			
 			*overlaps = *overlaps + 1;
		}
	}
}


/*------------------------------------------------------------------
 Find the particle pairs which overlap or violate the minimum
 gap specification
 -------------------------------------------------------------------*/ 
void find_repulsive_force(double xcm[npart][3],int overlap_pairs[npart][npart],
    double min_gap[npart][npart], double FR[npart][3])
{
	int ipart,jpart;
	int i,j;
	double xi[3],xj[3];
	double x,y,z,r;
	double num_min_gap=0.005; // numerical minimum gap
	double gap;
	double sign_x,sign_y;
	
	/* Initialize repulsive force to zero */
	for(ipart=0;ipart<npart;ipart++)
	{
		for(i=0;i<3;i++)
		{
			FR[ipart][i]=0.0;
		}
	}
	
	for(ipart=0;ipart<npart-1;ipart++)
	{
		for(jpart=ipart+1;jpart<npart;jpart++)
		{
			if (overlap_pairs[ipart][jpart] != 0)
			{
				for(i=0;i<3;i++)
				{
					xi[i] = xcm[ipart][i];
					xj[i] = xcm[jpart][i];
				}
		
				// difference vector
				x = xi[0]-xj[0];
				y = xi[1]-xj[1];
				z = xi[2]-xj[2];
				r = sqrt(x*x+y*y+z*z); 
				
				gap = min_gap[ipart][jpart];
				
				/* set a numerical minimum gap */
				if (gap < num_min_gap)
				{
					gap = num_min_gap;
				}
				
				/* correct for periodicity */
				if (fabs(x) > Lx/2) {
					sign_x = -1.0;
				}
				else {
					sign_x = 1.0;
				}
				
				if (fabs(y) > Ly/2) {
					sign_y = -1.0;
				}
				else {
					sign_y = 1.0;
				}
				
				FR[jpart][0] = (1.0 - GAP_MIN/gap)*x/r*sign_x;
				FR[jpart][1] = (1.0 - GAP_MIN/gap)*y/r*sign_y;
				FR[jpart][2] = (1.0 - GAP_MIN/gap)*z/r;
				
				FR[ipart][0] = -FR[jpart][0];
				FR[ipart][1] = -FR[jpart][1];
				FR[ipart][2] = -FR[jpart][2];
			}
		}
	}
	
}

/* -------------------------------------------------------------------
 * pushes the particle apart without rotating them to correct the 
 * overlaps 
 --------------------------------------------------------------------*/
void correct_overlaps(double xb[npart][beads][3], double nrm_b[npart][beads][3], double xcm[npart][3])
{
	int ipart,ibead,jpart;
	int i,j;
	int overlaps,overlaps0;
	double FR[npart][3];
	int overlap_pairs[npart][npart];
	double min_gap[npart][npart];
	double xcm_new[npart][3];
 	double REX=0.001;
	int itr_count=0;	
	FILE *fp;

	do 
	{	

		/* find the minimum gap violations */
		find_overlaps_gap(xb,nrm_b,overlap_pairs,min_gap, &overlaps);
		
		if (itr_count == 0)
		{
			overlaps0 = overlaps;
			/* set overlap indicator */
			for(ipart=0;ipart<npart;ipart++)
			{
				overlap_indicator[ipart]=0;
				for(jpart=0;jpart<npart;jpart++)
				{
					overlap_indicator[ipart] += overlap_pairs[ipart][jpart];
				}
			}
		}
		
		if (overlaps == 0 || itr_count > 1000)
		{
			
			break;
		}
		else
		{
#if (PRINT==1)
			printf("Number of overlaps = %d\t initial = %d\t iteration count = %d\n ",overlaps,overlaps0,itr_count);
#endif
			
		}
		
		/* find the repsulsive force */
		find_repulsive_force(xcm,overlap_pairs,min_gap,FR);
		
		/* move the particles */
		for(ipart=0;ipart<npart;ipart++)
		{
			for(i=0;i<3;i++)
			{
				xcm_new[ipart][i] = xcm[ipart][i] + REX*FR[ipart][i];
			}
		}
		
		/* find new xb */
		for(ipart=0;ipart<npart;ipart++)
		{
			for(ibead=0;ibead<beads;ibead++)
			{
				for(i=0;i<3;i++)
				{
					xb[ipart][ibead][i] = xb[ipart][ibead][i] + (xcm_new[ipart][i] - xcm[ipart][i]);
				}
			}
		}
		
		/* update center of mass */
		for(ipart=0;ipart<npart;ipart++)
		{
			for(i=0;i<3;i++)
			{
				xcm[ipart][i] = xcm_new[ipart][i];
			}
		}
		
		itr_count++;
		
	} while (1);
	
	fp = fopen("overlap.txt","a+");
	fprintf(fp,"%d\t%d\t%d\n ",overlaps,overlaps0,itr_count);
	fclose(fp);
}

/*-----------------------------------------------------------------------------------

-------------------------------------------------------------------------------------*/
void compute_total_grid_vel(double xb[npart][beads][3], double fb[npart][beads][3], int psteps)
{

  int i,j,k;
  int ibead,itriangle,jtriangle;
  double  xi[3];
  double r1;
  int i1,i2,i3;
  double pi;
  double u[3];
  int member;
  double weights[npart][beads];
  double uxf[nx][ny][nz], uyf[nx][ny][nz], uzf[nx][ny][nz],pf[nx][ny][nz];
  double ubc1[nx][ny][3], ubc2[nx][ny][3];
  double gaussx[nx][ny][nz],gaussy[nx][ny][nz],gaussz[nx][ny][nz],gaussz_p[nx][ny][nz];
  FILE *fp;
  char name[20];
  int ix,iy,iz;
  double dx,dy,dz;
  int ipart,jpart;
	double dudx[nx][ny][nz],dudy[nx][ny][nz],dudz[nx][ny][nz];
	double dvdx[nx][ny][nz],dvdy[nx][ny][nz],dvdz[nx][ny][nz];
	double dwdx[nx][ny][nz],dwdy[nx][ny][nz],dwdz[nx][ny][nz];

  pi = 4.0*atan(1.0);

  /* Mesh size */
  dx = Lx/nx;
  dy = Ly/ny;
  dz = Lz/(nz-1.0);



#if (OMP == 1)
omp_set_num_threads(NTHREADS);
#endif




  /*---------------- Global solution calculation -----------------*/

  /*----- Boundary Conditions ---- */
  /* Initialize */
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      for(k=0;k<3;k++)
      {
        ubc1[i][j][k]=0.0;
        ubc2[i][j][k]=0.0;
      }
    }
  }


  /* set ubc = - u_l for poiseuille flow */
  /* Top plate */
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = Lz;

      for(itriangle=0;itriangle<triangle_count_t[ix][iy];itriangle++)
      {
        jpart = triangle_list_t[ix][iy][itriangle][0];
        jtriangle = triangle_list_t[ix][iy][itriangle][1];

        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc2[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }

  /* bottom plate */
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      xi[0] = ix*dx;
      xi[1] = iy*dy;
      xi[2] = 0.0;

      for(itriangle=0;itriangle<triangle_count_b[ix][iy];itriangle++)
      {
        jpart = triangle_list_b[ix][iy][itriangle][0];
        jtriangle = triangle_list_b[ix][iy][itriangle][1];

        i1 = triangle[jtriangle][0];
        i2 = triangle[jtriangle][1];
        i3 = triangle[jtriangle][2];
        member=0; // performing non-singular integral as of now

        if (member == 1)
        {
          integrate2(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }
        else
        {
          integrate1(xb[jpart][i1],xb[jpart][i2],xb[jpart][i3],xi,fb[jpart][i1],fb[jpart][i2],fb[jpart][i3],u,alpha);
        }


        for(i=0;i<3;i++)
        {
          ubc1[ix][iy][i] -= u[i]/8/pi;
        }

      }
    }
  }

  /* Find weights for the global solution */
  /* Initialize weights for the global solution */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<beads;i++)
    {
      weights[ipart][i]=0.0;
    }
  }
  for(ipart=0;ipart<npart;ipart++)
  {
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];

      /* Store weights for the global calculation */
      weights[ipart][i1] += area[ipart][itriangle]/3.0;
      weights[ipart][i2] += area[ipart][itriangle]/3.0;
      weights[ipart][i3] += area[ipart][itriangle]/3.0;

    }
  }

  /* Distribute density to mesh */
  /* Initialize mesh force density */
  for(ix=0;ix<nx;ix++)
  {
    for(iy=0;iy<ny;iy++)
    {
      for(iz=0;iz<nz;iz++)
      {
        gaussx[ix][iy][iz] = 0.0;
        gaussy[ix][iy][iz] = 0.0;
        gaussz[ix][iy][iz] = 0.0; 
				gaussz_p[ix][iy][iz] = 0.0;
      }
    }
  }

  for(ipart=0;ipart<npart;ipart++)
  {
    distribute_density(xb[ipart],fb[ipart],gaussx,gaussy,gaussz,gaussz_p,weights[ipart],alpha,xz);
  }

	/* solve for the velocity and pressure at the mesh points */
	/* solve for the velocity and pressure at the mesh points */
	global_velocity_inhomogeneous(gaussx,gaussy,gaussz,gaussz_p,ubc1,ubc2,uxf,uyf,uzf,pf,
																dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,															
																xxi,cz,
																u1H,v1H,w1H,dwdz1H,dudz1H,dvdz1H,p1H,
															  u2H,v2H,w2H,dwdz2H,dudz2H,dvdz2H,p2H);



 /* Add local contribution to the velocity at grid points: doing in an umoptimized
  fashion, as this is for test purposes only */
 
 #if (OMP == 1)
 #pragma omp parallel for private(ix,iy,iz,xi,itriangle,ipart,i1,i2,i3,member,u)
 #endif
 for(ix=0;ix<nx;ix++)
 {
   for(iy=ny/2;iy<ny/2+1;iy++) // midplane only
   {
     for(iz=0;iz<nz;iz++)
     {
       xi[0] = ix*dx;
       xi[1] = iy*dy;
       xi[2] = iz*dz;
       
       for(ipart=0;ipart<npart;ipart++)
       {
         for(itriangle=0;itriangle<ntriangles;itriangle++)
         {                    
           i1 = triangle[itriangle][0];
           i2 = triangle[itriangle][1];
           i3 = triangle[itriangle][2];
           member=0; // performing non-singular integral as of now
           
           if (member == 1)
           {
             integrate2(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xi,fb[ipart][i1],fb[ipart][i2],fb[ipart][i3],u,alpha);
           }
           else
           {
             integrate1(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xi,fb[ipart][i1],fb[ipart][i2],fb[ipart][i3],u,alpha);
           }

           /* sum to grid velocity array */
           uxf[ix][iy][iz] -= u[0]/8/pi;
           uyf[ix][iy][iz] -= u[1]/8/pi;
           uzf[ix][iy][iz] -= u[2]/8/pi;
           
         }
       }
     }
   }
 }
 

  if (CONDOR==0)
  {
    sprintf(name,"output/vel_grid_%03d.vtk",psteps);
  }
  else
  {
    sprintf(name,"vel_grid_%03d.vtk",psteps);
  }

  fp = fopen(name,"w+");

  for(ix=0;ix<nx;ix++)
  {
    for(iz=0;iz<nz;iz++)
    {
      fprintf(fp,"%E\t%E\t%E\t%E\t%E\t%E\n",ix*dx,ny/2*dy,iz*dz,uxf[ix][ny/2][iz],uyf[ix][ny/2][iz],uzf[ix][ny/2][iz]);
    }
  }
  fclose(fp);



}


/*------------------------------------------------------------
Computes initial volume
------------------------------------------------------------*/
void compute_vol(double xb[beads][3], double xcm[3], double nrm[beads][3], double *Vol0, int ipart)
{
  
  int itriangle;
  int i1,i2,i3;
  
  
  /*------- Volume averaged velocity inside the capsule----------- */
  for(itriangle=0;itriangle<ntriangles;itriangle++)
  {
    i1 = triangle[itriangle][0];
    i2 = triangle[itriangle][1];
    i3 = triangle[itriangle][2];
    
		// normal of beads
    //integrate_vol(xb[i1],xb[i2],xb[i3],xcm,nrm[i1],nrm[i2],nrm[i3],Vol0);
		
		// normal of triangel
		integrate_vol(xb[i1],xb[i2],xb[i3],xcm,normal_t[ipart][itriangle],normal_t[ipart][itriangle],
									normal_t[ipart][itriangle],Vol0);
  }
  
}

/*----------------------------------------------*/
/* -------- randomly rotates the swimmers ------ */
/*----------------------------------------------*/
void rand_rotate(double xb[npart][beads][3], double xcm0[npart][3])
{
	int i,j,k,ipart;
	double Qrot[3][3],urot[3],norm_rot,xrot[3],xrotn[3];


  /* apply a random rotation */
  srand(time(NULL));
  
  
	for(ipart=0;ipart<npart;ipart++)
  {
		urot[0]= 1.0;//rand();
		urot[1]= 0.0;//rand();
		urot[2]= 1.0;//rand();
		norm_rot = sqrt(urot[0]*urot[0]+urot[1]*urot[1]+urot[2]*urot[2]);  
  
		urot[0] = urot[0]/norm_rot;
		urot[1] = urot[1]/norm_rot;
		urot[2] = urot[2]/norm_rot;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				if (i==j)
				{
					Qrot[i][j] = 1.0-2.0*urot[i]*urot[j];
				}
				else
				{
					Qrot[i][j]=-2.0*urot[i]*urot[j];
				}
			}
		}
 
    for(i=0;i<beads;i++)
    {
      
      for(j=0;j<3;j++)
      {
        xrot[j] = xb[ipart][i][j] - xcm0[ipart][j];
      }
      
      for(j=0;j<3;j++)
      {
        xrotn[j]=0.0;
        for(k=0;k<3;k++)
        {
          xrotn[j] += Qrot[j][k]*xrot[k];
        }
      }
      
      for(j=0;j<3;j++)
      {
        xb[ipart][i][j] = xcm0[ipart][j] + xrotn[j];
      }
    }
  }

}


/*------------------------------------------------------- */
/* -------- rotates the swimmers by directed angle------ */
/*--------------------------------------------------------*/
void dir_rotate(double xb[npart][beads][3], double xcm0[npart][3], double theta_rot)
{
	int i,j,k,ipart;
	double Qrot[3][3],theta,pi,xrot[3],xrotn[3];

	pi = 4.0*atan(1.0);
  /* apply a random rotation */
  // srand(time(NULL));
    theta = pi*(theta_rot/180);
    //printf("theta = %E\t cos = %E\t THETA_rot = %E\n",theta,cos(theta),theta_rot);
	for(ipart=0;ipart<npart;ipart++)
  {
		//urot[0]= 1.0;//rand();
		//urot[1]= 0.0;//rand();
		//urot[2]= 1.0;//rand();
		//norm_rot = sqrt(urot[0]*urot[0]+urot[1]*urot[1]+urot[2]*urot[2]);  
  
		//urot[0] = urot[0]/norm_rot;
		//urot[1] = urot[1]/norm_rot;
		//urot[2] = urot[2]/norm_rot;
		//for(i=0;i<3;i++)
		//{
			//for(j=0;j<3;j++)
			//{
				//if (i==j)
				//{
					//Qrot[i][j] = 1.0-2.0*urot[i]*urot[j];
				//}
				//else
				//{
					//Qrot[i][j]=-2.0*urot[i]*urot[j];
				//}
			//}
		//}
		
		/* First rotation by 90 about z-axis and theta about x-axis */
		//Qrot[0][0] = 1;
		//Qrot[0][1] = 0;
		//Qrot[0][2] = 0;
		//Qrot[1][0] = 0;
		//Qrot[1][1] = cos(theta);
		//Qrot[1][2] = -sin(theta);
		//Qrot[2][0] = 0;
		//Qrot[2][1] = sin(theta);
		//Qrot[2][2] = cos(theta);
		
		/* First rotation by 90 about z-axis and theta about y-axis(commented part) */
		Qrot[0][0] = cos(theta);
		Qrot[0][1] = 0;
		Qrot[0][2] = sin(theta);
		Qrot[1][0] = 0;
		Qrot[1][1] = 1;
		Qrot[1][2] = 0;
		Qrot[2][0] = -sin(theta);
		Qrot[2][1] = 0;
		Qrot[2][2] = cos(theta);
 
    for(i=0;i<beads;i++)
    {
      
      for(j=0;j<3;j++)
      {
        xrot[j] = xb[ipart][i][j] - xcm0[ipart][j];
      }
      
      for(j=0;j<3;j++)
      {
        xrotn[j]=0.0;
        for(k=0;k<3;k++)
        {
          xrotn[j] += Qrot[j][k]*xrot[k];
        }
      }
      
      for(j=0;j<3;j++)
      {
        xb[ipart][i][j] = xcm0[ipart][j] + xrotn[j];
      }
    }
  }

}


/*-----------------------------------------------------------------------------------
// rotates so as to make ibead=0  on the -ve z-axis
-------------------------------------------------------------------------------------*/
void rotate_z_align(double xb[beads][3],double xcm[3])
{
  int i,j,k;
  int ibead;
	double norm;
  double u[3],Q[3][3],zp[3];
	double xb_rot[beads][3];
	
	/* find the displacement wrt to the center of mass */
	for(ibead=0;ibead<beads;ibead++)
	{
		for(i=0;i<3;i++)
		{
			xb_rot[ibead][i] =  xb[ibead][i] - xcm[i];
		}
	}
		
  /* Store ibead=0 as the z-axis in the old frame of reference */
  zp[0] = xb_rot[0][0];
  zp[1] = xb_rot[0][1];
  zp[2] = xb_rot[0][2];

  /* find the rotation matrix: Use householder's algorithm */
  if (fabs(zp[2]-1) > 1E-10)
  {
		u[0] = 0.0 - zp[0];
    u[1] = 0.0 - zp[1];
    u[2] = 1.0 - zp[2];

    norm = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        
    /* normalize */
    u[0] = u[0]/norm;
    u[1] = u[1]/norm;
    u[2] = u[2]/norm;

    /* Find the Q matrix */
    for(i=0;i<3;i++)
    {
			for(j=0;j<3;j++)
			{
				if (i==j)
        {
					Q[i][j] = 1.0 - 2.0*u[i]*u[j];
        }
        else
        {
					Q[i][j] = -2.0*u[i]*u[j];
        }
      }
    }
   }
   else // no rotation necessary Q=I
   {
			for(i=0;i<3;i++)
      {
				for(j=0;j<3;j++)
        {
					if (i==j)
          {
						Q[i][j] = 1.0;
          }
          else
          {
						Q[i][j] = 0.0;
          }
        }
      }
    }

    /* transform the coordinates */
		for(ibead=0;ibead<beads;ibead++)
		{
			for(i=0;i<3;i++)
			{
				xb[ibead][i] = xcm[i];
				
				for(j=0;j<3;j++)
				{
					xb[ibead][i] += Q[i][j]*xb_rot[ibead][j];
				}
			}
		}
		
  
}

/*---------------------------------------------------------------*/
/*---- Computes the particle contribution to the stress tensor ---*/
/*---------------------------------------------------------------*/
void compute_bulk_stress_tensor(double xb[npart][beads][3],double xcm[npart][3],
     double fm[npart][beads][3],double ub[npart][beads][3],double nrm[npart][beads][3],
     double Sigma[npart][3][3],double Sigma_g[3][3],double Sigma_g1[3][3],
     double Sigma_g2[3][3],double Vcell, double viscr[npart], double Stime)
{
  int i,j,k,ibead,ipart;
  int itriangle;
  int i1,i2,i3;
  double S12,N1,N2;
  double S12a,N1a,N2a;
  double S12b,N1b,N2b;
  int npart1,npart2;
  FILE *fp;
  char name[20];
  double ndensity;

  npart2 = NRATIO*npart;
  npart1 = npart - npart2;

  
  /* Initialize */
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      Sigma_g[i][j]=0.0;
      Sigma_g1[i][j]=0.0;
      Sigma_g2[i][j]=0.0;
    }
  }
  
  /*------- Particle's contribution to bulk stress tensor----------- */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        Sigma[ipart][i][j]=0.0;
      }
    }
    
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];
			
			// normal of beads
   //   integrate_stress(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xcm[ipart],ub[ipart][i1],ub[ipart][i2],
   //                       ub[ipart][i3],nrm[ipart][i1],nrm[ipart][i2],nrm[ipart][i3],fm[ipart][i1],
   //                       fm[ipart][i2],fm[ipart][i3],Sigma[ipart],viscr[ipart]);
			
			// normal of triangle
            integrate_stress(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xcm[ipart],ub[ipart][i1],ub[ipart][i2],
                    ub[ipart][i3],normal_t[ipart][itriangle],normal_t[ipart][itriangle],normal_t[ipart][itriangle],
                    fm[ipart][i1],
                    fm[ipart][i2],fm[ipart][i3],Sigma[ipart],viscr[ipart]);
			
    }
  }


  /* find the particle's contribution to bulk stress */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        Sigma_g[i][j] += Sigma[ipart][i][j];
      }
    }
  }

  /* find the particle's contribution to bulk stress: Ist type */
  for(ipart=0;ipart<npart1;ipart++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        Sigma_g1[i][j] += Sigma[ipart][i][j];
      }
    }
  }

  /* find the particle's contribution to bulk stress: 2nd type */
  for(ipart=npart1;ipart<npart;ipart++)
  {
    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
      {
        Sigma_g2[i][j] += Sigma[ipart][i][j];
      }
    }
  }

/* Volume averaged contribution */
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      Sigma_g[i][j]=Sigma_g[i][j]/Vcell;
      Sigma_g1[i][j]=Sigma_g1[i][j]/Vcell*npart/npart1;
      if (npart2 > 0)
      {
        Sigma_g2[i][j]=Sigma_g2[i][j]/Vcell*npart/npart2;
      }
    }
  }

  S12 = 0.5*(Sigma_g[0][2]+Sigma_g[2][0]); // xy
  N1  = Sigma_g[0][0]-Sigma_g[2][2];   // xx - yy
  N2  = Sigma_g[2][2]-Sigma_g[1][1];   // yy - zz


  S12a = 0.5*(Sigma_g1[0][2]+Sigma_g1[2][0]);
  N1a  = Sigma_g1[0][0]-Sigma_g1[2][2];
  N2a  = Sigma_g1[2][2]-Sigma_g1[1][1];

  if (npart2 > 0)
  {
    S12b = 0.5*(Sigma_g2[0][2]+Sigma_g2[2][0]);
    N1b  = Sigma_g2[0][0]-Sigma_g2[2][2];
    N2b  = Sigma_g2[2][2]-Sigma_g2[1][1];
  }
  

  /* print */
  if (CONDOR == 1)
  {
    fp = fopen("Stress.txt","a+");
  }
  else
  {
    fp = fopen("output/Stress.txt","a+");
  }

  if(npart2 > 0)
  {
    fprintf(fp,"%E %E %E\t%E %E %E\t%E %E %E\n",S12,N1,N2,S12a,N1a,N2a,S12b,N1b,N2b);
  }
  else
  {
    fprintf(fp,"%E %E %E\t%E %E %E\n",S12,N1,N2,S12a,N1a,N2a);
  }

  
  fclose(fp);

  if (PRINT == 1)
  {
    printf("S12 = %E\t N1 = %E\t N2 = %E\n",S12,N1,N2);
    printf("S12a = %E\t N1a = %E\t N2a = %E\n",S12a,N1a,N2a);
    if (npart2>0)
    {
      printf("S12b = %E\t N1b = %E\t N2b= %E\n",S12b,N1b,N2b);
    }
  }

  /* print per particle stress tensor */
  ndensity = npart/Vcell;
  for(ipart=0;ipart<npart;ipart++)
  {
    if (CONDOR==1)
    {
      sprintf(name,"stress_%03d.txt",ipart);
    }
    else
    {
      sprintf(name,"output/stress_%03d.txt",ipart);
    }
    S12 = 0.5*(Sigma[ipart][0][2]+Sigma[ipart][2][0]); // xy
    N1  = Sigma[ipart][0][0]-Sigma[ipart][2][2];   // xx - yy
    N2  = Sigma[ipart][2][2]-Sigma[ipart][1][1];   // yy - zz
    
    fp = fopen(name,"a+");
    fprintf(fp,"%E\t%E\t%E\n",S12*ndensity,N1*ndensity,N2*ndensity);
    fclose(fp);
  }
    
  
}
/*-------------------------------------------------------------------------------------*/
/* ----- Computes volume averaged velocity and rate of rotation of particles ---------*/
/*-------------------------------------------------------------------------------------*/
void volume_avg_u_omega(double ub[npart][beads][3],double nrm[npart][beads][3],
                        double xcm[npart][3],double xb[npart][beads][3],double Stime)
{

  int i,j,k,ibead,ipart;
  int itriangle;
  int i1,i2,i3;
  double pi;
  double uvol[npart][3],uvolw[npart];
  double xcm_v[npart][3];
  double Omega[3][3],omega[npart][3];
  int epsilon[3][3][3];
  double x1[3],x2[3],x3[3];
  FILE *fp;
  char name[20];
  
  pi = 4.0*atan(1.0);


  /* set permutation operator */
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      for(k=0;k<3;k++)
      {
        if(i==j || i==k || j==k)
        {
          epsilon[i][j][k]=0;
        }
        else
        {
          if(i==0)
          {
            if(j==1)
            {
              epsilon[i][j][k]=1;
            }
            else
            {
              epsilon[i][j][k]=-1;
            }
          }
          if(i==1)
          {
            if(j==2)
            {
              epsilon[i][j][k]=1;
            }
            else
            {
              epsilon[i][j][k]=-1;
            }
          }
          if(i==2)
          {
            if(j==0)
            {
              epsilon[i][j][k]=1;
            }
            else
            {
              epsilon[i][j][k]=-1;
            }
          }
        }
      }
    }
  }
  
  /*------- Volume averaged velocity and rotation inside the capsule----------- */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<3;i++)
    {
      uvol[ipart][i]=0.0;
      for(j=0;j<3;j++)
      {
        Omega[i][j]=0.0;
      }
    }
    
    uvolw[ipart]=0.0;
    
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];
			
			// normal of beads
      //integrate_trans_vol(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xcm[ipart],ub[ipart][i1],ub[ipart][i2],ub[ipart][i3],nrm[ipart][i1],nrm[ipart][i2],nrm[ipart][i3],uvol[ipart],Omega,&uvolw[ipart]);
			
			// normal of triangle
			integrate_trans_vol(xb[ipart][i1],xb[ipart][i2],xb[ipart][i3],xcm[ipart],ub[ipart][i1],ub[ipart][i2],ub[ipart][i3],normal_t[ipart][itriangle],normal_t[ipart][itriangle],normal_t[ipart][itriangle],uvol[ipart],Omega,&uvolw[ipart]);
    }

    /* volume averaged velocity */
    for(i=0;i<3;i++)
    {
      uvol[ipart][i] = uvol[ipart][i]/uvolw[ipart];
    }

    /* volume averaged rotation rate */
    for(i=0;i<3;i++)
    {
      omega[ipart][i]=0.0;
      for(j=0;j<3;j++)
      {
        for(k=0;k<3;k++)
        {
          omega[ipart][i] += 0.5*epsilon[i][j][k]*Omega[j][k];
        }
      }
    }
  }


  /*------- Volume averaged center of mass of the capsule----------- */
  for(ipart=0;ipart<npart;ipart++)
  {
    for(i=0;i<3;i++)
    {
      xcm_v[ipart][i]=0.0;
    }
    for(itriangle=0;itriangle<ntriangles;itriangle++)
    {
      i1 = triangle[itriangle][0];
      i2 = triangle[itriangle][1];
      i3 = triangle[itriangle][2];

      for(i=0;i<3;i++)
      {
        x1[i] = xb[ipart][i1][i] - xcm[ipart][i];
        x2[i] = xb[ipart][i2][i] - xcm[ipart][i];
        x3[i] = xb[ipart][i3][i] - xcm[ipart][i];
      }

			// normal at beads
      //integrate_rsq_vol(x1,x2,x3,nrm[ipart][i1],nrm[ipart][i2],nrm[ipart][i3],xcm_v[ipart]);
      
      // normal of triangle
      integrate_rsq_vol(x1,x2,x3,normal_t[ipart][itriangle],
			normal_t[ipart][itriangle],normal_t[ipart][itriangle],xcm_v[ipart]);

    }
    for(i=0;i<3;i++)
    {
      xcm_v[ipart][i] = xcm_v[ipart][i]/uvolw[ipart] + xcm[ipart][i];
    }
  }





  if(PRINT==1)
  {
    if (npart==1)
    {
      printf("Volume average U: %E\t%E\t%E\n",uvol[0][0],uvol[0][1],uvol[0][2]);
      printf("Vol avg CM: %E\t%E\t%E\n",xcm_v[0][0],xcm_v[0][1],xcm_v[0][2]);
      printf("Volume = %E\n",uvolw[0]/4/pi*3);
    }
  }

  for(ipart=0;ipart<npart;ipart++)
  {
    if (CONDOR==1)
    {
      sprintf(name,"pos_vel_%03d.txt",ipart);
    }
    else
    {
      sprintf(name,"output/pos_vel_%03d.txt",ipart);
    }
      
    fp = fopen(name,"a+");
    fprintf(fp,"%E\t",Stime);
#if (DETAILED==1)    
    fprintf(fp,"%E\t%E\t%E\t",xcm_v[ipart][0],xcm_v[ipart][1],xcm_v[ipart][2]);
    fprintf(fp,"%E\t%E\t%E\t",uvol[ipart][0],uvol[ipart][1],uvol[ipart][2]);
    fprintf(fp,"%E\t%E\t%E\n",omega[ipart][0],omega[ipart][1],omega[ipart][2]);
#else
		//fprintf(fp,"%E\t%E\t%E\n",xcm_v[ipart][0],xcm_v[ipart][1],xcm_v[ipart][2]);
		 fprintf(fp,"%E\t%E\t%E\n",uvol[ipart][0],uvol[ipart][1],uvol[ipart][2]);
#endif
		

    fclose(fp);
  }

}

/*-----------------------------------------------------------------------------------
Find Taylor deformation parameter
-------------------------------------------------------------------------------------*/
void deformation_parameter(double xb[beads][3],double xcm[3],double Le[3], double thetae[2],double phie[2])
{
  int i,j;
  double max=-100,min=100;
  double dist;
  double x,y,z;


  for(i=0;i<beads;i++)
  {
    x = xb[i][0]-xcm[0];
    y = xb[i][1]-xcm[1];
    z = xb[i][2]-xcm[2];
    dist = sqrt(x*x+y*y+z*z);

    if (dist < min)
    {
      min = dist;
      Le[0]    = dist;
      thetae[0] = acos(y/dist);
      phie[0]   = atan(z/x);
    }

    if (dist > max)
    {
      max = dist;
      Le[1]    = dist;
      thetae[1] = acos(y/dist);
      phie[1]   = atan(z/x);
    }
  }


}

/*-----------------------------------------------------------------------------------
Find Taylor deformation parameter using mass moment of inertia tensor
-------------------------------------------------------------------------------------*/
void mmoi_deformation_parameter(double xb[beads][3], double xcm0[3], double DD[2], double Stime, int ipart, int print)
{
	int i, j, k, ii, jj, kk;
	double mmoi[3][3],delta[3][3], xcm[3];
	double dxij[3], dxjk[3], dxik[3];
	double xcm_to_cm[3], normal[3], v[3];
	double rij,rjk,rik,s,area_l,norm,theta_cm_to_cm, r2;
	int N, LDA, LWORK, INFO;
	double eigenv[3];
	double* work;
	double wkopt;
	double lmax, lmin, lmid, pi;
	FILE *fp;
	char name[20];

	pi = 4*atan2(1,1);
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			mmoi[i][j] = 0.0;
			delta[i][j] = 0.0;
		}
	}
	
	for(i=0;i<3;i++)
	{
		delta[i][i] = 1.0;
	}
	
	for(i=0;i<ntriangles;i++)
	{
		for(j=0;j<3;j++)
		{
			dxij[j]  = xb[triangle[i][1]][j]-xb[triangle[i][0]][j];
			dxjk[j]  = xb[triangle[i][2]][j]-xb[triangle[i][1]][j];
			dxik[j]  = xb[triangle[i][0]][j]-xb[triangle[i][2]][j];
		}

		rij = sqrt(dot_product(3,dxij,dxij));
		rjk = sqrt(dot_product(3,dxjk,dxjk));
		rik = sqrt(dot_product(3,dxik,dxik));

		s = 0.5*(rij+rjk+rik);

		area_l = sqrt(s*(s-rij)*(s-rjk)*(s-rik));

		/* center of mass of the triangle */
		for(j=0;j<3;j++)
		{
			xcm[j] = (xb[triangle[i][0]][j] + xb[triangle[i][1]][j] + xb[triangle[i][2]][j])/3.0;
		}
	
		cross_product(dxij, dxjk, normal);
    
		norm = sqrt(dot_product(3,normal,normal));
		for(j=0;j<3;j++)
		{
			normal[j] = normal[j]/norm;
		}
    
		for(j=0;j<3;j++)
		{
			xcm_to_cm[j] = xcm0[j] - xcm[j];
		}
    
    
		norm = sqrt(dot_product(3,xcm_to_cm,xcm_to_cm));
		for(j=0;j<3;j++)
		{
			xcm_to_cm[j] = xcm_to_cm[j]/norm;
		}
    
    
		theta_cm_to_cm = acos(dot_product(3,normal,xcm_to_cm));
    
		if(theta_cm_to_cm < pi/2.0)
		{
			for(j=0;j<3;j++)
			{
				normal[j] = -normal[j];
			}
		}
	  
		for(j=0;j<3;j++)
		{
			v[j] = xcm[j] - xcm0[j];
		}
		r2 = dot_product(3,v,v);
		
	    for(ii=0;ii<3;ii++)
	    {
			for(jj=0;jj<3;jj++)
			{
				for(kk=0;kk<3;kk++)
				{
					mmoi[ii][jj] += (r2*v[kk]*delta[ii][jj] - v[ii]*v[jj]*v[kk])*normal[kk]*area_l;
				}
			}
		}
	}
	
	mmoi[1][0] = 0.0; mmoi[2][0] = 0.0; mmoi[2][1] = 0.0;
  	N=3;LDA=3;LWORK=-1;
    dsyev_( "V", "L", &N, mmoi, &LDA, eigenv, &wkopt, &LWORK, &INFO );
    LWORK = (int)wkopt;
    printf("%d",LWORK);
    work = (double*)malloc( LWORK*sizeof(double) );
    /* Solve eigenproblem */
    dsyev_( "V", "L", &N, mmoi, &LDA, eigenv, work, &LWORK, &INFO );
    /* Check for convergence */
     if( INFO > 0 ) {
          printf( "The algorithm failed to compute eigenvalues.\n" );
          exit( 1 );
       }
  	 	  	
	//kaiser(mmoi, 3, 3, eigenv, trace, sume, ier); 
 
	lmax = sqrt(eigenv[2]+eigenv[1]-eigenv[0]);
	lmin = sqrt(eigenv[1]+eigenv[0]-eigenv[2]);
	lmid = sqrt(eigenv[0]+eigenv[2]-eigenv[1]);
  
	DD[0] = (lmax-lmin)/(lmax+lmin);
	DD[1] = atan(mmoi[0][2]/mmoi[0][0]);
	
	/* Print Eigenvalues and Eigenvectors */
    if (PRINT == 1 || print == 1)
    {
      if (CONDOR==1)
        {
          sprintf(name,"Eigenval_vector.txt");
        }
        else
        {
          sprintf(name,"output/Eigenval_vector.txt");
        }

        fp=fopen(name,"a+");
        fprintf(fp,"%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Stime,ipart,eigenv[0],eigenv[1],eigenv[2],
				mmoi[0][0],mmoi[0][1],mmoi[0][2],mmoi[1][0],mmoi[1][1],mmoi[1][2],mmoi[2][0],mmoi[2][1],mmoi[2][2]);
		fclose(fp);
	}
  
}

/*-----------------------------------------------------------------------------------
 find area weighted center of mass of the body
-------------------------------------------------------------------------------------*/
void compute_center_mass_area(double xb[beads][3], double xcm0[3])
{
  int i,j;
  double xcm[3],xcm_to_cm[3];
  double dxij[3],dxjk[3],dxik[3];
  double rij,rjk,rik,s,theta_cm_to_cm;
  double area_sum=0.0;
  double area_l;
  
  for(i=0;i<3;i++)
  {
    xcm0[i]=0.0;
  }
  for(i=0;i<ntriangles;i++)
  {
    for(j=0;j<3;j++)
    {
      dxij[j]  = xb[triangle[i][1]][j]-xb[triangle[i][0]][j];
      dxjk[j]  = xb[triangle[i][2]][j]-xb[triangle[i][1]][j];
      dxik[j]  = xb[triangle[i][0]][j]-xb[triangle[i][2]][j];
    }

    rij = sqrt(dot_product(3,dxij,dxij));
    rjk = sqrt(dot_product(3,dxjk,dxjk));
    rik = sqrt(dot_product(3,dxik,dxik));

    s = 0.5*(rij+rjk+rik);

    area_l = sqrt(s*(s-rij)*(s-rjk)*(s-rik));
    area_sum += area_l;

    /* center of mass of the triangle */
    for(j=0;j<3;j++)
    {
      xcm[j] = (xb[triangle[i][0]][j] + xb[triangle[i][1]][j] + xb[triangle[i][2]][j])/3.0;
    }

    /* sum to center of mass of the body */
    for(j=0;j<3;j++)
    {
      xcm0[j] += xcm[j]*area_l;
    }
    
  }

  /* Find CM */
  for(j=0;j<3;j++)
  {
    xcm0[j] = xcm0[j]/area_sum;
  }

  if(PRINT==1)
  {
    if (npart == 1)
    {
      printf("xcm=%lf\tycm=%lf\tzcm=%lf\trad=%lf\n",xcm0[0],xcm0[1],xcm0[2],sqrt(area_sum/4/3.14));
    }
  }

}

/*-----------------------------------------------------------------------------------
computes the traction at nodes due to the imposed flow:finf
-------------------------------------------------------------------------------------*/
void compute_finf(double finf[beads][3],double xb[beads][3],double nrm[beads][3])
{
  int ibead;
  double pi;
  double U0;
  double gdot=GDOT;

  /* Compute the stress tensor and then traction at each node */
  if (SHEAR==0)
  {
    /* set U0=centerline velocity */
    U0 = gdot*Lz/4.0;
    for(ibead=0;ibead<beads;ibead++)
    {
      finf[ibead][0] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][0] + 4*U0/Lz*(1-2*xb[ibead][2]/Lz)*nrm[ibead][2];
      finf[ibead][1] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][1];
      finf[ibead][2] = 8*U0*xb[ibead][0]/Lz/Lz*nrm[ibead][2] + 4*U0/Lz*(1-2*xb[ibead][2]/Lz)*nrm[ibead][0];
    }
  }
  else
  {
    for(ibead=0;ibead<beads;ibead++)
    {
      finf[ibead][0] = gdot*nrm[ibead][2];
      finf[ibead][1] = 0.0;
      finf[ibead][2] = gdot*nrm[ibead][0];
    }
  }
}

/*--------------------------------------------------*/
/* sets the viscosity ratio and membrane stiffness */
/*-------------------------------------------------*/
void set_viscr(double viscr[npart], double krbc_r[npart])
{
    int ipart;
    double randr;
    int npart1,npart2;
    
    npart2 = NRATIO*npart;
    npart1 = npart - npart2;
    /* set viscosity ratio */
    if (LAMBDA == 1) {
        for(ipart=0;ipart<npart;ipart++) {
            viscr[ipart]=1.0;
        }
    } 
    else {
        for(ipart=0;ipart<npart1;ipart++) {
      	    viscr[ipart]=LAMBDA1;
        }
        for(ipart=npart1;ipart<npart;ipart++) {
        		viscr[ipart]=LAMBDA2;
        }
    }
    /* set membrane shear modulus  */
    for(ipart=0;ipart<npart1;ipart++) {
      	krbc_r[ipart]   = 1.0;
    }
    for(ipart=npart1;ipart<npart;ipart++) {
      	krbc_r[ipart] = KRBC_R;
    }
}

