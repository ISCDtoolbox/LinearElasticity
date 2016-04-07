#include "elastic.h"
#include "sparse.h"


/* find boundary conds in list */
pCl getCl(pSol sol,int ref,int elt) {
  pCl     pcl;
  int     i;

  for (i=0; i<sol->nbcl; i++) {
    pcl = &sol->cl[i];
    if ( (pcl->ref == ref) && (pcl->elt == elt) )  return(pcl);
  }
  return(0);
}

/* retrieve physical properties in list */
int getMat(pSol sol,int ref,double *lambda,double *mu) {
  pMat   pm;
  int    i;

  *lambda = LS_LAMBDA;
  *mu     = LS_MU;
  if ( sol->nmat == 0 )  return(1);
  for (i=0; i<sol->nmat; i++) {
    pm = &sol->mat[i];
    if ( pm->ref == ref ) {
      *lambda = pm->lambda;
      *mu     = pm->mu;
      return(1);
    }
  }

  return(0);
}

/* triangle area */
static inline double area_2d(double *a,double *b,double *c) {
  return(0.5 * ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])));
}

/* return length of edge + outer normal to segment */ 
static double length(double *a,double *b,double n[2]) {
  double   dd;

  n[0] = -(b[1] - a[1]);
  n[1] =   b[0] - a[0];
  dd   = sqrt(n[0]*n[0] + n[1]*n[1]);
  if ( dd > LS_EPSD ) {
    n[0] /= dd;
    n[1] /= dd;
  }
  return(dd);
}

static int setTGV_2d(LSst *lsst,pCsr A) {
	pCl      pcl;
  pPoint   ppt;
  pEdge    pa;
  int      k,i,dof;

  /* at vertices */
  if ( lsst->sol.clelt & LS_ver ) {
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( pcl && pcl->typ == Dirichlet ) {
        csrSet(A,2*(k-1)+0,2*(k-1)+0,LS_TGV);
        csrSet(A,2*(k-1)+1,2*(k-1)+1,LS_TGV);
      }
    }
  }
  /* at edge nodes */
  else if ( lsst->sol.clelt & LS_edg ) {
    dof = lsst->info.typ == P1 ? 2 : 3;
    for (k=1; k<=lsst->info.na; k++) {
      pa  = &lsst->mesh.edge[k];
      pcl = getCl(&lsst->sol,pa->ref,LS_edg);
      if ( pcl && pcl->typ == Dirichlet ) {
        for (i=0; i<dof; i++) {
          csrSet(A,2*(pa->v[i]-1)+0,2*(pa->v[i]-1)+0,LS_TGV);
          csrSet(A,2*(pa->v[i]-1)+1,2*(pa->v[i]-1)+1,LS_TGV);
        }
      }
    }
  }

  return(1);
}


static pCsr matA_P1_2d(LSst *lsst) {
  pCsr     A;
  pTria    pt;
  double   Ae[6][6],DeD[4][4],m[2][2],mm[4][6],nn[4][6],*a,*b,*c;
  double   lambda,mu,det,idet,area;
  int      nr,nc,nbe,i,j,k,s,ia,ja,ig,jg,il,ic;
  
  /* memory allocation (rough estimate) */
  nr  = nc = 2*lsst->info.np;
  nbe = 10*lsst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  memset(mm,0,4*6*sizeof(double));
  memset(DeD,0,16*sizeof(double));

  /* store values in A */
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];

    /* tD E D */
    if ( !getMat(&lsst->sol,pt->ref,&lambda,&mu) )  continue;
    DeD[0][0] = DeD[3][3] = 2.0*mu + lambda;
    DeD[0][3] = DeD[3][0] = lambda;
    DeD[1][1] = DeD[1][2] = DeD[2][1] = DeD[2][2] = mu; 

    /* measure of K */
    a = &lsst->mesh.point[pt->v[0]].c[0];
    b = &lsst->mesh.point[pt->v[1]].c[0];
    c = &lsst->mesh.point[pt->v[2]].c[0];

    /* m = tBT^-1 */
    det  = (b[1]-c[1])*(a[0]-c[0]) - (a[1]-c[1])*(b[0]-c[0]);
    idet = 1.0 / det;
    m[0][0] = idet*(b[1]-c[1]);    m[0][1] = idet*(c[1]-a[1]);
    m[1][0] = idet*(c[0]-b[0]);    m[1][1] = idet*(a[0]-c[0]);

    /* mm = (tBT^-1) Dp */
    mm[0][0] = mm[2][3] = m[0][0];
    mm[0][1] = mm[2][4] = m[0][1];
    mm[0][2] = mm[2][5] = -(m[0][0]+m[0][1]);
    mm[1][0] = mm[3][3] = m[1][0];
    mm[1][1] = mm[3][4] = m[1][1];
    mm[1][2] = mm[3][5] = -(m[1][0]+m[1][1]);

    /* nn = DeD mm */
    for (i=0; i<4; i++) {
      for (j=0; j<6; j++) {
        nn[i][j] = 0.0;
        for (s=0; s<4; s++)
          nn[i][j] += DeD[i][s] * mm[s][j];
      }
    }

    area = 0.5 * det;
    memset(Ae,0,6*6*sizeof(double));
    /* Ae = tmm * nn */
    for (i=0; i<6; i++) {
      for (j=i; j<6; j++) {
        for (s=0; s<4; s++) {
          Ae[i][j] += area * mm[s][i] * nn[s][j];
        }
      }
    }

    /* stifness matrix */
    for (i=0; i<6; i++) {
      ig = pt->v[i % 3];
      ia = 2*(ig-1) + (i / 3);
      for (j=i; j<6; j++) {
        if ( fabs(Ae[i][j]) < LS_EPSD )  continue;
        jg = pt->v[j % 3];
        ja = 2*(jg-1) + (j / 3);
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
  setTGV_2d(lsst,A);
  csrPack(A);

  if ( lsst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}


static pCsr matA_P2_2d(LSst *lsst) {
  return(0);
}


/* build right hand side vector and set boundary conds. */
static double *rhsF_P1_2d(LSst *lsst) {
  pTria    pt;
  pEdge    pa;
  pPoint   ppt;
  pCl      pcl;
  double  *F,*vp,area,lon,n[2],w[2],*a,*b,*c;
  int      i,k,nc;

  if ( lsst->info.verb == '+' )  fprintf(stdout,"     gravity and body forces\n");
  F = (double*)calloc(lsst->info.dim * lsst->info.np,sizeof(double));
  assert(F);

  /* gravity as external force */
  if ( lsst->sol.cltyp & Gravity ) {
    nc = 0;
    for (k=1; k<=lsst->info.nt; k++) {
      pt = &lsst->mesh.tria[k];

      /* measure of K */
      a = &lsst->mesh.point[pt->v[0]].c[0]; 
      b = &lsst->mesh.point[pt->v[1]].c[0]; 
      c = &lsst->mesh.point[pt->v[2]].c[0]; 
      area = area_2d(a,b,c) / 3.0;
      for (i=0; i<3; i++) {
        F[2*(pt->v[i]-1)+0] += area * lsst->sol.gr[0];
        F[2*(pt->v[i]-1)+1] += area * lsst->sol.gr[1];
      }
      nc++;
    }
    if ( lsst->info.verb == '+' )  fprintf(stdout,"     %d gravity values assigned\n",nc);
  }

  /* nodal boundary conditions */
  if ( lsst->sol.clelt & LS_ver ) {
    nc = 0;
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        vp = pcl->att == 'f' ? &lsst->sol.u[2*(k-1)] : &pcl->u[0];
        F[2*(k-1)+0] = LS_TGV * vp[0];
        F[2*(k-1)+1] = LS_TGV * vp[1];
      }
      else if ( pcl->typ == Load ) {
        vp = pcl->att == 'f' ? &lsst->sol.u[2*(k-1)] : &pcl->u[0];
        F[2*(k-1)+0] = vp[0];
        F[2*(k-1)+1] = vp[1];
      }
      nc++;
    }
    if ( lsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d nodal values\n",nc);
  }

  /* external load along boundary edges */
  if ( lsst->sol.clelt & LS_edg ) {
    nc = 0;
    for (k=1; k<=lsst->info.na; k++) {
      pa  = &lsst->mesh.edge[k];
      pcl = getCl(&lsst->sol,pa->ref,LS_edg);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        for (i=0; i<2; i++) {
          vp = pcl->att == 'f' ? &lsst->sol.u[2*(pa->v[i]-1)] : &pcl->u[0];
          F[2*(pa->v[i]-1)+0] = LS_TGV * vp[0];
          F[2*(pa->v[i]-1)+1] = LS_TGV * vp[1];
        }
        nc++;
      }
      /* load along normal direction (normal component) */
      else if ( pcl->typ == Load ) {
        a = &lsst->mesh.point[pa->v[0]].c[0];
        b = &lsst->mesh.point[pa->v[1]].c[0];
        lon = length(a,b,n);
        if ( pcl->att == 'n' ) {
          w[0] = 0.5 * lon * pcl->u[0] * n[0];
          w[1] = 0.5 * lon * pcl->u[0] * n[1];
        }
        else {
          vp = pcl->att == 'f' ? &lsst->sol.u[2*(k-1)] : &pcl->u[0];
          w[0] = 0.5 * lon * vp[0];
          w[1] = 0.5 * lon * vp[1];
        }
        F[2*(pa->v[0]-1)+0] += w[0];
        F[2*(pa->v[0]-1)+1] += w[1];
        F[2*(pa->v[1]-1)+0] += w[0];
        F[2*(pa->v[1]-1)+1] += w[1];
        nc++;
      }
    }
    if ( lsst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d load values\n",nc);
  }

  return(F);
}


/* build right hand side vector and set boundary conds. */
static double *rhsF_P2_2d(LSst *lsst) {
  return(0);
}


/* 2d linear elasticity */
int elasti1_2d(LSst *lsst) {
  pCsr     A;
  int      ier;
  char     stim[32];

  /* -- Part I: matrix assembly */
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");

  /* counting P2 nodes (for dylib) */
	if ( lsst->info.typ == P2 && !lsst->info.np2 ) {
		lsst->info.np2 = hashar_2d(lsst);
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

  /* build matrix */
  A = lsst->info.typ == P1 ? matA_P1_2d(lsst) : matA_P2_2d(lsst);
  lsst->sol.F = lsst->info.typ == P1 ? rhsF_P1_2d(lsst) : rhsF_P2_2d(lsst);

  /* free mesh structure + boundary conditions */
  if ( !lsst->info.xport && lsst->info.mfree ) {
		free(lsst->mesh.tria);
    if ( !lsst->info.zip )  free(lsst->mesh.point);
	}

  /* -- Part II: solver */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Solving linear system:");  fflush(stdout);
    ier = csrPrecondGrad(A,lsst->sol.u,lsst->sol.F,&lsst->sol.res,&lsst->sol.nit,1);
    if ( ier <= 0 )
      fprintf(stdout,"\n # convergence problem: %d\n",ier);
    else
      fprintf(stdout," %E in %d iterations\n",lsst->sol.res,lsst->sol.nit);
	}
  else {
    ier = csrPrecondGrad(A,lsst->sol.u,lsst->sol.F,&lsst->sol.res,&lsst->sol.nit,1);
  }

  /* free memory */
  csrFree(A);
  free(lsst->sol.F);

  return(ier > 0);
}

