#include "elastic.h"
#include "ls_calls.h"
#include "sparse.h"


/* compute triangle area and unit normal in 3d */
static double area_3d(double *a,double *b,double *c,double *n) {
  double    ux,uy,uz,vx,vy,vz,dd,dd1;

  ux = b[0] - a[0];
  uy = b[1] - a[1];
  uz = b[2] - a[2];

  vx = c[0] - a[0];
  vy = c[1] - a[1];
  vz = c[2] - a[2];

  n[0] = uy*vz - uz*vy;
  n[1] = uz*vx - ux*vz;
  n[2] = ux*vy - uy*vx;
  dd   = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( dd > LS_EPSD ) {
    dd1 = 1.0 / dd;
    n[0] *= dd1;
    n[1] *= dd1;
    n[2] *= dd1;
  }

  return(0.5*dd);
}


/* compute volume of tetra */
static inline double volume(double *a,double *b,double *c,double *d) {
  double  ax,ay,az,bx,by,bz,vol;

  ax = b[0] - a[0];
  ay = b[1] - a[1];
  az = b[2] - a[2];
  bx = c[0] - a[0];
  by = c[1] - a[1];
  bz = c[2] - a[2];
  vol = (d[0]-a[0]) * (ay*bz - az*by) + (d[1]-a[1]) * (az*bx - ax*bz) \
      + (d[2]-a[2]) * (ax*by - ay*bx);

  return(fabs(vol) / 6.0);
}


/* invert 3x3 non-symmetric matrix */
static int invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < LS_EPSD )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}


static int setTGV_3d(LSst *lsst,pCsr A) {
  pCl      pcl;
  pTria    ptt;
  pPoint   ppt;
  int      k,ig,nbpt;
  char     i;

  /* at vertices */
  if ( lsst->sol.clelt & LS_ver ) {
    for (k=1; k<=lsst->info.np+lsst->info.np2; k++) {
      ppt = &lsst->mesh.point[k];
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( pcl && pcl->typ == Dirichlet ) {
        csrSet(A,3*(k-1)+0,3*(k-1)+0,LS_TGV);
        csrSet(A,3*(k-1)+1,3*(k-1)+1,LS_TGV);
        csrSet(A,3*(k-1)+2,3*(k-1)+2,LS_TGV);
      }
    }
  }
  /* at elements nodes (triangles) */
  if ( lsst->sol.clelt & LS_tri ) {
		nbpt = lsst->info.typ == P1 ? 3 : 6;
		for (k=1; k<=lsst->info.nt; k++) {
      ptt = &lsst->mesh.tria[k];
      if ( !ptt->v[0] )  continue;
      pcl = getCl(&lsst->sol,ptt->ref,LS_tri);
      if ( pcl && pcl->typ == Dirichlet ) {
        for (i=0; i<nbpt; i++) {
					csrSet(A,3*(ptt->v[i]-1)+0,3*(ptt->v[i]-1)+0,LS_TGV);
					csrSet(A,3*(ptt->v[i]-1)+1,3*(ptt->v[i]-1)+1,LS_TGV);
					csrSet(A,3*(ptt->v[i]-1)+2,3*(ptt->v[i]-1)+2,LS_TGV);
        }
      }
    } 
  }

  return(1);
}


/* build stifness matrix P2 */
int assMat_P2(pTetra pt,pCsr A,double *DeD,double *im,double vol) {
  double  Ae[30][30],mm[9][30],nn[9][30],Dp[3][10],wp;
  int     i,j,p,s,ia,ja,ic,il,ig,jg;

  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                         {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

  memset(Ae,0,30*30*sizeof(double));

  /* boucle sur les 5 points de quadrature */
  for (p=0; p<5; p++) {
    /* Dp */
    Dp[0][0]=4*q[p][0]-1; Dp[0][1]=0;            Dp[0][2]=0;           Dp[0][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[1][0]=0;           Dp[1][1]=4*q[p][1]-1;  Dp[1][2]=0;           Dp[1][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[2][0]=0;           Dp[2][1]=0;            Dp[2][2]=4*q[p][2]-1; Dp[2][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    
    Dp[0][4]=4*q[p][1];   Dp[0][5]=4*q[p][2];    Dp[0][6]=4*(1-2*q[p][0]-q[p][1]-q[p][2]); 
    Dp[1][4]=4*q[p][0];   Dp[1][5]=0;            Dp[1][6]=-4*q[p][0]; 
    Dp[2][4]=0;           Dp[2][5]=4*q[p][0];    Dp[2][6]=-4*q[p][0];   
  
    Dp[0][7]=0;           Dp[0][8]=-4*q[p][1];                         Dp[0][9]=-4*q[p][2];
    Dp[1][7]=4*q[p][2];   Dp[1][8]=4*(1-q[p][0]-2*q[p][1]-q[p][2]);    Dp[1][9]=-4*q[p][2];
    Dp[2][7]=4*q[p][1];   Dp[2][8]=-4*q[p][1];                         Dp[2][9]=4*(1-q[p][0]-q[p][1]-2*q[p][2]);
    
    /* mm = (tBt^-1) Dp */
    memset(mm,0,9*30*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<10; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
        mm[i+3][j+10] = mm[i][j];
        mm[i+6][j+20] = mm[i][j];
      }
    }
   
    /* nn = DeD mm */
    memset(nn,0,9*30*sizeof(double));
    for (i=0; i<9; i++) {
      for (j=0; j<30; j++) {
        for (s=0; s<9; s++)
          nn[i][j] += DeD[i*9+s] * mm[s][j];
      }
    }
  
    /* Ae = vol tmm DeD mm */
    wp = vol * w[p];
    for (i=0; i<30; i++) {
      for (j=i; j<30; j++) {
        for (s=0; s<9; s++)
          Ae[i][j] += wp * mm[s][i] * nn[s][j];
      }
    }
  }
  
  /* stifness matrix */
  for (i=0; i<30; i++) {
    ig = pt->v[i % 10];
    ia = 3*(ig-1) + (i / 10);
    for (j=i; j<30; j++) {
      if ( fabs(Ae[i][j]) < LS_EPSD )  continue;
      jg = pt->v[j % 10];
      ja = 3*(jg-1) + (j / 10);
      il = LS_MIN(ia,ja);
      ic = LS_MAX(ia,ja);
      for (s=A->row[il]; s<A->row[il+1]; s++) {
        if ( A->col[s] == -1 ) {
          A->col[s] = ic;
          A->val[s] = Ae[i][j];
          break;
        }
        else if ( A->col[s] == ic ) {
          A->val[s] += Ae[i][j];
          break;
        }
      }
      assert(s < A->row[il+1]);
    }
  }

  return(1);
}

static pCsr matA_P2_3d(LSst *lsst) {
  pCsr     A;
  pTetra   pt;
  double  *a,*b,*c,*d,DeD[81],m[9],im[9],Ae[30][30],mm[9][30],nn[9][30],Dp[3][10];
  double   wp,lambda,mu,vol;
  int      i,j,k,p,s,ia,ja,il,ic,ig,jg,nr,nc,nbe;
  static double w[5]    = { -4./5., 9./20., 9./20., 9./20.,  9./20. };
  static double q[5][3] = { {1./4.,1./4.,1./4.}, {1./2.,1./6.,1./6.}, {1./6.,1./2.,1./6.}, 
                            {1./6.,1./6.,1./2.}, {1./6.,1./6.,1./6.} };

	/* memory allocation (rough estimate) */
	nr  = nc = 3 * (lsst->info.np+lsst->info.np2);
  nbe = 30 * lsst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  memset(DeD,0,81*sizeof(double));

  /* Dp for quadrature points */
  for (p=0; p<5; p++) {
    /* Dp */
    Dp[0][0]=4*q[p][0]-1;  Dp[0][1]=0;            Dp[0][2]=0;            Dp[0][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[1][0]=0;            Dp[1][1]=4*q[p][1]-1;  Dp[1][2]=0;            Dp[1][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    Dp[2][0]=0;            Dp[2][1]=0;            Dp[2][2]=4*q[p][2]-1;  Dp[2][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
                           
    Dp[0][4]=4*q[p][1];    Dp[0][5]=4*q[p][2];    Dp[0][6]=4*(1-2*q[p][0]-q[p][1]-q[p][2]); 
    Dp[1][4]=4*q[p][0];    Dp[1][5]=0;            Dp[1][6]=-4*q[p][0]; 
    Dp[2][4]=0;            Dp[2][5]=4*q[p][0];    Dp[2][6]=-4*q[p][0];   
                           
    Dp[0][7]=0;            Dp[0][8]=-4*q[p][1];                       Dp[0][9]=-4*q[p][2];
    Dp[1][7]=4*q[p][2];    Dp[1][8]=4*(1-q[p][0]-2*q[p][1]-q[p][2]);  Dp[1][9]=-4*q[p][2];
    Dp[2][7]=4*q[p][1];    Dp[2][8]=-4*q[p][1];                       Dp[2][9]=4*(1-q[p][0]-q[p][1]-2*q[p][2]);
  }

  /* Fill stiffness matrix A */
  for (k=1; k<=lsst->info.ne; k++) {
    pt = &lsst->mesh.tetra[k];
    if ( !pt->v[0] )  continue;

    /* tD E D (related to material properties) */
    if ( !getMat(&lsst->sol,pt->ref,&lambda,&mu) )  continue;
    DeD[0]  = DeD[40] = DeD[80] = 2.0 * mu + lambda;
    DeD[4]  = DeD[8]  = DeD[36] = DeD[44] = DeD[72] = DeD[76] = lambda;
    DeD[10] = DeD[12] = DeD[20] = DeD[24] = DeD[28] = DeD[30] = mu; 
    DeD[50] = DeD[52] = DeD[56] = DeD[60] = DeD[68] = DeD[70] = mu;

    /* measure of K */
    a = &lsst->mesh.point[pt->v[0]].c[0]; 
    b = &lsst->mesh.point[pt->v[1]].c[0]; 
    c = &lsst->mesh.point[pt->v[2]].c[0]; 
    d = &lsst->mesh.point[pt->v[3]].c[0]; 

    /* mm = tB^-1 */
    for (i=0; i<3; i++) {
      m[i+0] = a[i] - d[i]; // jacobian
      m[i+3] = b[i] - d[i];
      m[i+6] = c[i] - d[i];
    }
    if ( !invmatg(m,im) )  return(0);
    vol = volume(a,b,c,d);

    /* mm = (tBt^-1) Dp */
    /* discrete strain element e_mm  e_ij = d_ju_i + d_iu_j */
	  memset(Ae,0,30*30*sizeof(double));

    /* mm = (tBt^-1) Dp */
    memset(mm,0,9*30*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<10; j++) {
        for (s=0; s<3; s++) 
          mm[i][j]   += im[i*3+s] * Dp[s][j];
        mm[i+3][j+10] = mm[i][j];
        mm[i+6][j+20] = mm[i][j];
      }
    }

    /* nn = DeD mm */
    memset(nn,0,9*30*sizeof(double));
    for (i=0; i<9; i++) {
      for (j=0; j<30; j++) {
        for (s=0; s<9; s++)
          nn[i][j] += DeD[i*9+s] * mm[s][j];
      }
    }

	  /* loop over 5 quadrature points */
	  for (p=0; p<5; p++) {
	    /* Dp */
	    Dp[0][0]=4*q[p][0]-1; Dp[0][1]=0;            Dp[0][2]=0;           Dp[0][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
	    Dp[1][0]=0;           Dp[1][1]=4*q[p][1]-1;  Dp[1][2]=0;           Dp[1][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
	    Dp[2][0]=0;           Dp[2][1]=0;            Dp[2][2]=4*q[p][2]-1; Dp[2][3]=4*(q[p][0]+q[p][1]+q[p][2])-3; 
    
	    Dp[0][4]=4*q[p][1];   Dp[0][5]=4*q[p][2];    Dp[0][6]=4*(1-2*q[p][0]-q[p][1]-q[p][2]); 
	    Dp[1][4]=4*q[p][0];   Dp[1][5]=0;            Dp[1][6]=-4*q[p][0]; 
	    Dp[2][4]=0;           Dp[2][5]=4*q[p][0];    Dp[2][6]=-4*q[p][0];   
  
	    Dp[0][7]=0;           Dp[0][8]=-4*q[p][1];                         Dp[0][9]=-4*q[p][2];
	    Dp[1][7]=4*q[p][2];   Dp[1][8]=4*(1-q[p][0]-2*q[p][1]-q[p][2]);    Dp[1][9]=-4*q[p][2];
	    Dp[2][7]=4*q[p][1];   Dp[2][8]=-4*q[p][1];                         Dp[2][9]=4*(1-q[p][0]-q[p][1]-2*q[p][2]);
    
	    /* Ae = vol tmm DeD mm */
	    wp = vol * w[p];
	    for (i=0; i<30; i++) {
	      for (j=i; j<30; j++) {
	        for (s=0; s<9; s++)
	          Ae[i][j] += wp * mm[s][i] * nn[s][j];
	      }
	    }
	  }

	  /* stifness matrix */
	  for (i=0; i<30; i++) {
	    ig = pt->v[i % 10];
	    ia = 3*(ig-1) + (i / 10);
	    for (j=i; j<30; j++) {
	      if ( fabs(Ae[i][j]) < LS_EPSD )  continue;
	      jg = pt->v[j % 10];
	      ja = 3*(jg-1) + (j / 10);
        if ( ia < ja ) {
          il = ia;
          ic = ja;
        }
        else {
          il = ja;
          ic = ia;
        }
				/* a(il,jl)= Ae[i][j]) */
        csrPut(A,il,ic,Ae[i][j]);
	    }
	  }
	}

  /* Set large value for Dirichlet conditions */
  setTGV_3d(lsst,A);
	csrPack(A);

  return(A);
}


static pCsr matA_P1_3d(LSst *lsst) {
  pCsr     A;
  pTetra   pt;
  double  *a,*b,*c,*d,DeD[81],m[9],im[9],Ae[12][12],mm[9][12],nn[9][12],Dp[3][4];
  double   lambda,mu,vol;
  int      i,j,k,s,ia,ja,il,ic,ig,jg,nr,nc,nbe;

	/* memory allocation (rough estimate) */
	nr  = nc = 3 * lsst->info.np;
  nbe = 12 * lsst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  memset(DeD,0,81*sizeof(double));

  /* Dp for quadrature points */
  Dp[0][0]=1;  Dp[0][1]=0;  Dp[0][2]=0;  Dp[0][3]=-1; 
  Dp[1][0]=0;  Dp[1][1]=1;  Dp[1][2]=0;  Dp[1][3]=-1; 
  Dp[2][0]=0;  Dp[2][1]=0;  Dp[2][2]=1;  Dp[2][3]=-1; 

  /* Fill stiffness matrix A */
  for (k=1; k<=lsst->info.ne; k++) {
    pt = &lsst->mesh.tetra[k];

    /* tD E D */
    if ( !getMat(&lsst->sol,pt->ref,&lambda,&mu) )  continue;
    DeD[0]  = DeD[40] = DeD[80] = 2.0 * mu + lambda;
    DeD[4]  = DeD[8]  = DeD[36] = DeD[44] = DeD[72] = DeD[76] = lambda;
    DeD[10] = DeD[12] = DeD[20] = DeD[24] = DeD[28] = DeD[30] = mu; 
    DeD[50] = DeD[52] = DeD[56] = DeD[60] = DeD[68] = DeD[70] = mu;

    /* measure of K */
    a = &lsst->mesh.point[pt->v[0]].c[0]; 
    b = &lsst->mesh.point[pt->v[1]].c[0]; 
    c = &lsst->mesh.point[pt->v[2]].c[0]; 
    d = &lsst->mesh.point[pt->v[3]].c[0]; 

    /* mm = tB^-1 (Jacobian) */
    for (i=0; i<3; i++) {
      m[i+0] = a[i] - d[i];
      m[i+3] = b[i] - d[i];
      m[i+6] = c[i] - d[i];
    }
    if ( !invmatg(m,im) )  return(0);
    vol = volume(a,b,c,d);

    /* mm = (tBt^-1) Dp */
    /* discrete strain element e_mm  e_ij = d_ju_i + d_iu_j */
    memset(mm,0,9*12*sizeof(double));
    for (i=0; i<3; i++) {
      for (j=0; j<4; j++) {
        for (s=0; s<3; s++)
          mm[i][j]   += im[i*3+s] * Dp[s][j]; /* inv(J)*Dp from the derivative */
        mm[i+3][j+4] = mm[i][j];
        mm[i+6][j+8] = mm[i][j];
      }
    }

    /* nn = DeD mm */
    /* discrete stress element */
    for (i=0; i<9; i++) {
      for (j=0; j<12; j++) {
        nn[i][j] = 0.0;
        for (s=0; s<9; s++)
          nn[i][j] += DeD[i*9+s] * mm[s][j];
      }
    }

    /* Ae = vol tmm nn */
    /* stiffness matrix element K = e^t (De) */
    memset(Ae,0,12*12*sizeof(double));
    for (i=0; i<12; i++) {
      for (j=i; j<12; j++) {
        for (s=0; s<9; s++)
          Ae[i][j] += vol * mm[s][i] * nn[s][j];
      }
    }

    /* stifness matrix */
    for (i=0; i<12; i++) {
      ig = pt->v[i % 4];
      ia = 3*(ig-1) + (i / 4);
      for (j=i; j<12; j++) {
        if ( fabs(Ae[i][j]) < LS_EPSD )  continue;
        jg = pt->v[j % 4];
        ja = 3*(jg-1) + (j / 4);
        if ( ia < ja ) {
          il = ia;
          ic = ja;
        }
        else {
          il = ja;
          ic = ia;
        }
        csrPut(A,il,ic,Ae[i][j]);
      }
    }
  }
  setTGV_3d(lsst,A);
  csrPack(A);

  if ( lsst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}


/* build right hand side vector and set boundary conds. */
static double *rhsF_3d(LSst *lsst) {
  pTetra   pt;
	pTria    ptt;
  pPoint   ppt;
  pCl      pcl;
  double  *F,*vp,aire,vol,n[3],w[3],*a,*b,*c,*d;
  int      k,ig,nbpt,nc,size;
  char     i;
  
  if ( lsst->info.verb == '+' )  fprintf(stdout,"     gravity and body forces\n");
  size = lsst->info.dim * (lsst->info.np+lsst->info.np2);
  F = (double*)calloc(size,sizeof(double));
  assert(F);

  /* gravity as external force */
  if ( lsst->sol.cltyp & Gravity ) {
    nc = 0;
    for (k=1; k<=lsst->info.ne; k++) {
      pt = &lsst->mesh.tetra[k];

      /* measure of K */
      a = &lsst->mesh.point[pt->v[0]].c[0]; 
      b = &lsst->mesh.point[pt->v[1]].c[0]; 
      c = &lsst->mesh.point[pt->v[2]].c[0];
      d = &lsst->mesh.point[pt->v[3]].c[0];
      vol = volume(a,b,c,d) / 4.0;
      for (i=0; i<4; i++) {
        F[3*(pt->v[i]-1)+0] += vol * lsst->sol.gr[0];
        F[3*(pt->v[i]-1)+1] += vol * lsst->sol.gr[1];
        F[3*(pt->v[i]-1)+2] += vol * lsst->sol.gr[2];
      }
      nc++;
    }
    if ( lsst->info.verb == '+' )  fprintf(stdout,"     %d gravity values assigned\n",nc);
  }

  /* nodal boundary conditions */
  if ( lsst->sol.clelt & LS_ver ) {
    nc = 0;
    for (k=1; k<=lsst->info.np+lsst->info.np2; k++) {
      ppt = &lsst->mesh.point[k];
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        vp = pcl->att == 'f' ? &lsst->sol.u[3*(k-1)] : &pcl->u[0];
        F[3*(k-1)+0] = LS_TGV * vp[0];
        F[3*(k-1)+1] = LS_TGV * vp[1];
        F[3*(k-1)+2] = LS_TGV * vp[2];
      }
      else if ( pcl->typ == Load ) {
        vp = pcl->att == 'f' ? &lsst->sol.u[3*(k-1)] : &pcl->u[0];
        F[3*(k-1)+0] += vp[0];
        F[3*(k-1)+1] += vp[1];
        F[3*(k-1)+2] += vp[2];
      }
      nc++;
    }
    if ( lsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d nodal values\n",nc);
  }

  /* external load along boundary triangles */
  if ( lsst->sol.clelt & LS_tri ) {
    nbpt = lsst->info.typ == P1 ? 3 : 6;
    nc = 0;
    for (k=1; k<=lsst->info.nt; k++) {
      ptt = &lsst->mesh.tria[k];
      pcl = getCl(&lsst->sol,ptt->ref,LS_tri);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        for (i=0; i<nbpt; i++) {
          vp = pcl->att == 'f' ? &lsst->sol.u[3*(ptt->v[i]-1)] : &pcl->u[0];
          w[0] = LS_TGV * vp[0];
          w[1] = LS_TGV * vp[1];
          w[2] = LS_TGV * vp[2];
          F[3*(ptt->v[i]-1)+0] = w[0];
          F[3*(ptt->v[i]-1)+1] = w[1];
          F[3*(ptt->v[i]-1)+2] = w[2];
        }
        nc++;
      }
      else if ( pcl->typ == Load ) { 
				/* quadrature formula: to be modified eventually */
        a = &lsst->mesh.point[ptt->v[0]].c[0];
        b = &lsst->mesh.point[ptt->v[1]].c[0];
        c = &lsst->mesh.point[ptt->v[2]].c[0];
        aire = area_3d(a,b,c,n) / nbpt;
        if ( pcl->att == 'n' ) {
          w[0] = aire * pcl->u[0] * n[0];
          w[1] = aire * pcl->u[0] * n[1];
          w[2] = aire * pcl->u[0] * n[2];
        }
        else {
          vp = pcl->att == 'f' ? &lsst->sol.u[3*(k-1)] : &pcl->u[0];
          w[0] = aire * vp[0];
          w[1] = aire * vp[1];
          w[2] = aire * vp[2];
        }
        for (i=0; i<nbpt; i++) {
          F[3*(ptt->v[i]-1)+0] += w[0];
          F[3*(ptt->v[i]-1)+1] += w[1];
          F[3*(ptt->v[i]-1)+2] += w[2];
        }
        nc++;
      }
    }
    if ( lsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d load values\n",nc);
  }

	return(F);
}


/* 3d linear elasticity */
int elasti1_3d(LSst *lsst) {
  pCsr     A;
	int      ier;
	char     stim[32];
  const char typ[3] = {'0', '1', '2'};

  /* -- Part I: matrix assembly */
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");

  /* counting P2 nodes (for dylib) */
	if ( lsst->info.typ == P2 && !lsst->info.np2 ) {
		lsst->info.np2 = hashar_3d(lsst);
		if ( lsst->info.np2 == 0 ) {
			fprintf(stdout," # Error on P2 nodes.\n");
			return(0);
		}
	}

  /* allocating memory (for dylib) */
  if ( !lsst->sol.u ) {
    lsst->sol.u  = (double*)calloc(lsst->info.dim*(lsst->info.npi+lsst->info.np2),sizeof(double));
    assert(lsst->sol.u);
  }

  /* build matrix and right-hand side */
  A = lsst->info.typ == P1 ? matA_P1_3d(lsst) : matA_P2_3d(lsst);
  lsst->sol.F = lsst->info.typ == P1 ? rhsF_3d(lsst) : rhsF_3d(lsst);

  /* free mesh structure + boundary conditions */
  if ( !lsst->info.xport && lsst->info.mfree ) {
		free(lsst->mesh.tetra);
		if ( lsst->info.nt )    free(lsst->mesh.tria);
    if ( !lsst->info.zip )  free(lsst->mesh.point);
	}

  /* -- Part II: solver */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Solving linear system:");  fflush(stdout);
    ier = csrPrecondGrad(A,lsst->sol.u,lsst->sol.F,&lsst->sol.res,&lsst->sol.nit,0);
    if ( ier <= 0 )
      fprintf(stdout,"\n # convergence problem: %d\n",ier);
    else
      fprintf(stdout," %E in %d iterations\n",lsst->sol.res,lsst->sol.nit);
	}
  else {
    ier = csrPrecondGrad(A,lsst->sol.u,lsst->sol.F,&lsst->sol.res,&lsst->sol.nit,0);
  }

  /* free memory */
  csrFree(A);
  free(lsst->sol.F);

  return(ier > 0);
}

