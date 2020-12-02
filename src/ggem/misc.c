void init_normal_beads(double normal[NTRIANGLE_MAX][3], double nrm[beads][3], int ipart)
{
  double wx[beads],wy[beads],wz[beads];
  int i;
  int i1,i2,i3;
  int itriangle;
  double norm;
  
  /* Initialize weights for the global solution */
  for(i=0;i<beads;i++)
  {
    wx[i]=0.0;
    wy[i]=0.0;
    wz[i]=0.0;
  }
  
  
  
  /* Store weights for the global calculation */
  /* Also used for computing the area weighted normal average for each bead */
  for(itriangle=0;itriangle<ntriangles;itriangle++)
  {
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
  for(i=0;i<beads;i++)
  {
    norm = wx[i]*wx[i] + wy[i]*wy[i] + wz[i]*wz[i];
    norm = sqrt(norm);
    
    nrm[i][0] = wx[i]/norm;
    nrm[i][1] = wy[i]/norm;
    nrm[i][2] = wz[i]/norm;
  }
  

}


/*****************************************************************************/

void compute_normal_curvature(double nrm[beads][3], double curv[beads],
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
  }

  /* Remaining vertices with coordination number 6 */
  for(ibead=12;ibead<beads;ibead++) {
    /* store coordinates of the connected beads relative to ibead */
    for(i=0;i<6;i++) {
      jbead = beads_to_beads[ibead][i];
      for(j=0;j<3;j++) {
        xneigh[i][j] = xb[jbead][j]-xb[ibead][j];
      }
    }

    do {
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
        rdist2[i] = sqrt(x_tr[i][0]*x_tr[i][0] + x_tr[i][1]*x_tr[i][1]
                + x_tr[i][2]*x_tr[i][2]);
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
      for(i=0;i<3;i++) {
        nrm_g[i]=0.0;
        for(j=0;j<3;j++) {
          nrm_g[i] += Ct[j][i]*nrm_l[j];
        }
      }

      /* find error */
      error=0.0;
      for(i=0;i<3;i++) {
        error += pow(zp[i]-nrm_g[i],2);
      }
      error = sqrt(error);

      /* store in normal array */
      for(i=0;i<3;i++) {
        nrm[ibead][i] = nrm_g[i];
      }

    } while (error > 1E-3);

    /* store curvature */
    curv[ibead] = -brhs2[2] - brhs2[4];
  }
}

/*****************************************************************************/
/*****************************************************************************/
