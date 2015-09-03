#include "elastic.h"
#include "ls_calls.h"
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

  for (i=0; i<sol->nmat; i++) {
    pm = &sol->mat[i];
    if ( pm->ref == ref ) {
      *lambda = pm->lambda;
      *mu     = pm->mu;
      return(1);
    }
  }
  *lambda = LS_LAMBDA;
  *mu     = LS_MU;
  return(0);
}

/* triangle area */
static inline double area_2d(double *a,double *b,double *c) {
  double    ux,uy,vx,vy,dd;

  ux = b[0] - a[0];
  uy = b[1] - a[1];
  vx = c[0] - a[0];
  vy = c[1] - a[1];
  dd = 0.5 * (ux*vy - uy*vx);
  return(dd);
}

/* return length of edge + outer normal to segment */ 
static double length(double *a,double *b,double n[2]) {
  double   ax,ay,dd;

  ax = b[0] - a[0];
  ay = b[1] - a[1];
  n[0] = -ay;
  n[1] =  ax;
  dd   = sqrt(ax*ax + ay*ay);
  if ( dd > EPSD ) {
    n[0] *= 1.0 / dd;
    n[1] *= 1.0 / dd;
  }
  return(dd);
}

static int setTGV_2d(LSst *lsst,Hash *hash,pCsr A) {
	pCl      pcl;
  pPoint   ppt;
  int      k;

  /* at vertices */
  if ( lsst->sol.cltyp & LS_ver ) {
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( pcl && pcl->typ == Dirichlet ) {
        csrSet(A,2*(k-1)+0,2*(k-1)+0,LS_TGV);
        csrSet(A,2*(k-1)+1,2*(k-1)+1,LS_TGV);
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
    if ( !pt->v[0] )  continue;

    /* tD E D */
    if ( !getMat(&lsst->sol,pt->ref,&lambda,&mu) )  continue;
    DeD[0][0] = DeD[3][3] = 2*mu + lambda;
    DeD[0][3] = DeD[3][0] = lambda;
    DeD[1][1] = DeD[1][2] = DeD[2][1] = DeD[2][2] = mu; 

    /* measure of K */
    a = &lsst->mesh.point[pt->v[0]].c[0]; 
    b = &lsst->mesh.point[pt->v[1]].c[0]; 
    c = &lsst->mesh.point[pt->v[2]].c[0]; 

    /* m = tBT^-1 */
    det  = (b[1]-c[1])*(a[0]-c[0])-(a[1]-c[1])*(b[0]-c[0]);
    if ( det < EPSD )  continue;
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
        for (s=0; s<4; s++)
          Ae[i][j] += area * mm[s][i] * nn[s][j];
      }
    }

    /* stifness matrix */
    for (i=0; i<6; i++) {
      ig = pt->v[i % 3];
      ia = 2*(ig-1) + (i / 3);
      for (j=i; j<6; j++) {
        if ( fabs(Ae[i][j]) < EPSD )  continue;
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

  setTGV_2d(lsst,0,A);
  csrPack(A);
  if ( abs(lsst->info.imprim) > 5 || lsst->info.ddebug )
    fprintf(stdout,"     A: %6d x %6d  sparsity %7.4f%%\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}

/* build right hand side vector and set boundary conds. */
static double *rhsF_P1_2d(LSst *lsst) {
  pTria    pt;
  pEdge    pa;
  pPoint   ppt;
  pCl      pcl,pcl1;
  double  *F,*va,*vp,area,lon,n[2],w[2],*a,*b,*c;
  int      k,ig,size;
  char     i;

  if ( abs(lsst->info.imprim) > 5 )  fprintf(stdout,"     Gravity and body forces\n");
  size = lsst->info.dim * lsst->info.np;
  F = (double*)calloc(size,sizeof(double));
  assert(F);

  /* gravity as external force */
  if ( lsst->info.load && (1<<0) ) {
    for (k=1; k<=lsst->info.nt; k++) {
      pt = &lsst->mesh.tria[k];
      if ( !pt->v[0] )  continue;

      /* measure of K */
      a = &lsst->mesh.point[pt->v[0]].c[0]; 
      b = &lsst->mesh.point[pt->v[1]].c[0]; 
      c = &lsst->mesh.point[pt->v[2]].c[0]; 
      area = area_2d(a,b,c) / 3.0;
      for (i=0; i<3; i++) {
        ig = pt->v[i];
        F[2*(ig-1)+0] += area * lsst->info.gr[0];
        F[2*(ig-1)+1] += area * lsst->info.gr[1];
      }
    }
  }

  /* nodal boundary conditions */
  if ( lsst->sol.cltyp & LS_ver ) {
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      if ( !ppt->ref )  continue;
      pcl = getCl(&lsst->sol,ppt->ref,LS_ver);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        vp = pcl->att == 'f' ? &lsst->sol.bc[2*(k-1)] : &pcl->u[0];
        F[2*(k-1)+0] = LS_TGV * vp[0];
        F[2*(k-1)+1] = LS_TGV * vp[1];
      }
    }
  }

  /* external load along boundary edges (01/2011) */
  if ( lsst->sol.cltyp & LS_edg ) {
    for (k=1; k<=lsst->info.na; k++) {
      pa  = &lsst->mesh.edge[k];
      if ( !pa->ref )  continue;
      pcl = getCl(&lsst->sol,pa->ref,LS_edg);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        va = pcl->att == 'f' ? &lsst->sol.bc[2*(k-1)] : &pcl->u[0];
        for (i=0; i<2; i++) {
          ppt = &lsst->mesh.point[pa->v[i]];
          if ( pa->ref != ppt->ref ) {
            pcl1 = getCl(&lsst->sol,ppt->ref,LS_edg);
            vp = pcl1->att == 'f' ? &lsst->sol.bc[2*(pa->v[i]-1)] : &pcl1->u[0];
            F[2*(pa->v[i]-1)+0] = LS_TGV * vp[0];
            F[2*(pa->v[i]-1)+1] = LS_TGV * vp[1];
          }
          else {
            F[2*(pa->v[i]-1)+0] = LS_TGV * va[0];
            F[2*(pa->v[i]-1)+1] = LS_TGV * va[1];            
          }
        }
      }
      /* load along normal direction (normal component) */
      else if ( pcl->typ == Load ) { // formule de quadrature a changer eventuellement
        a = &lsst->mesh.point[pa->v[0]].c[0];
        b = &lsst->mesh.point[pa->v[1]].c[0];
        lon = length(a,b,n) / 2.0;
        w[0] = pcl->u[0] * n[0];
        w[1] = pcl->u[0] * n[1];
        for (i=0; i<2; i++) {
          ig  = pa->v[i];
          F[2*(ig-1)+0] += lon * w[0];
          F[2*(ig-1)+1] += lon * w[1];
        }
      }
    }
  }

  return(F);
}


/* 2d linear elasticity */
int elasti1_2d(LSst *lsst) {
  pCsr     A;
  double  *F;
  int      ier;
  char     stim[32];

  /* -- Part I: matrix assembly */
  chrono(ON,&lsst->info.ctim[3]);
  if ( abs(lsst->info.imprim) > 4 ) {
    fprintf(stdout,"  1.1 ASSEMBLY\n");
    fprintf(stdout,"     Assembly FE matrix\n");
  }

  /* counting P2 nodes (for dylib) */
	if ( lsst->info.typ == P2 && !lsst->info.np2 )  lsst->info.np2 = hashar(lsst);

  /* allocating memory (for dylib) */
  if ( !lsst->sol.u ) {
    lsst->sol.u  = (double*)calloc(lsst->info.dim * (lsst->info.npi+lsst->info.np2),sizeof(double));
    assert(lsst->sol.u);
  }

  /* build matrix */
  A = 0;
  F = 0;
	
	if ( lsst->info.typ == P1 ) {
    A = matA_P1_2d(lsst);
    F = rhsF_P1_2d(lsst);
	}
  else {
    /*
		A = matA_P2_2d(mesh,sol);
    F = rhsF_P2_2d(mesh,sol);
		*/
  }

  chrono(OFF,&lsst->info.ctim[3]);
  printim(lsst->info.ctim[3].gdif,stim);
  if ( abs(lsst->info.imprim) > 4 )  fprintf(stdout,"     [Time: %s]\n",stim);

  /* free mesh structure + boundary conditions */
  if ( lsst->info.mfree ) {
		free(lsst->mesh.tria);
    if ( !lsst->info.zip )  free(lsst->mesh.point);
	}
  if ( lsst->info.typ == P2 )  free(lsst->hash.item);

  /* -- Part II: solver */
  if ( abs(lsst->info.imprim) > 4 )  fprintf(stdout,"  1.2 SOLVING\n");
  chrono(ON,&lsst->info.ctim[4]);
  ier = csrPrecondGrad(A,lsst->sol.u,F,&lsst->sol.err,&lsst->sol.nit,1);
  chrono(OFF,&lsst->info.ctim[4]);
  if ( abs(lsst->info.imprim) > 0 ) {
    if ( ier <= 0 )
      fprintf(stdout,"  ## SOL NOT CONVERGED: ier= %d\n",ier);
    else if ( abs(lsst->info.imprim) > 4 ) {
      fprintf(stdout,"  %%%% CONVERGENCE: err= %E  nit= %d\n",lsst->sol.err,lsst->sol.nit);
		  printim(lsst->info.ctim[4].gdif,stim);
      fprintf(stdout,"     [Time: %s]\n",stim);
		}
	}

  /* free memory */
  csrFree(A);
  free(F);

  return(ier > 0);
}

