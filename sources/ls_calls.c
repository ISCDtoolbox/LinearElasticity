#include "elastic.h"
#include "ls_calls.h"


LSst *LS_init(int dim,int ver,char typ,char mfree) {
	LSst   *lsst;

  /* default values */
  lsst = (LSst*)calloc(1,sizeof(LSst));
  memset(&lsst->mesh,0,sizeof(Mesh));

  /* solution structure */
  memset(&lsst->sol,0,sizeof(Sol));
	lsst->sol.cl  = (Cl*)calloc(LS_CL,sizeof(Cl));
  lsst->sol.mat = (Mat*)calloc(LS_MAT,sizeof(Mat));
  lsst->sol.err = LS_RES;
  lsst->sol.nit = LS_MAXIT;

  /* global parameters */
  lsst->info.dim    = dim;
	lsst->info.ver    = ver;
  lsst->info.imprim = -99;
  lsst->info.ddebug = 0;
  lsst->info.zip    = 0;
  lsst->info.typ    = typ;
  lsst->info.mfree  = mfree;

  /* init timer */
  tminit(lsst->info.ctim,TIMEMAX);
  chrono(ON,&lsst->info.ctim[0]);

  return(lsst);
}


/* free global data structure */
int LS_stop(LSst *lsst) {
	char   stim[32];

	/* release memory */
  free(lsst->sol.bc);
	free(lsst->sol.cl);
	free(lsst->sol.mat);
  if ( lsst->sol.namein )  free(lsst->sol.namein);
	if ( lsst->sol.nameout ) free(lsst->sol.nameout);

  chrono(OFF,&lsst->info.ctim[0]);
  if ( abs(lsst->info.imprim) > 0 ) {
	  printim(lsst->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
  }

	return(1);
}

/* set params (facultative) */
void LS_setPar(LSst *lsst,int imp,int deb,int zip) {
	lsst->info.imprim = imp;
	lsst->info.ddebug = deb;
	lsst->info.zip    = zip;
}

/* handle boundary conditions:
  typ= P0, P1, P2
  ref= integer
  att= char 'v', 'f', 'n'
  elt= enum LS_ver, LS_edg, LS_tri, LS_tet */
int LS_setBC(LSst *lsst,int typ,int ref,char att,int elt,double *u) {
  int    i;

	if ( lsst->sol.nbcl == LS_CL-1 )  return(0);

  lsst->sol.nbcl++;
	lsst->sol.cl[lsst->sol.nbcl].typ = typ;
	lsst->sol.cl[lsst->sol.nbcl].ref = ref;
	lsst->sol.cl[lsst->sol.nbcl].att = att;
  lsst->sol.cl[lsst->sol.nbcl].elt = elt;

  if ( lsst->sol.cl[lsst->sol.nbcl].typ == Dirichlet ) {
	  if ( (lsst->sol.cl[lsst->sol.nbcl].att != 'v') && (lsst->sol.cl[lsst->sol.nbcl].att != 'f') ) {
      fprintf(stdout,"  %%%% Wrong format\n");
      return(0);
		}
	}
  else if ( lsst->sol.cl[lsst->sol.nbcl].typ == Load ) {
		if ( (lsst->sol.cl[lsst->sol.nbcl].att != 'v') && (lsst->sol.cl[lsst->sol.nbcl].att != 'f') \
		  && (lsst->sol.cl[lsst->sol.nbcl].att != 'n') ) {
      fprintf(stdout,"  %%%% Wrong format\n");
      return(0);
		}
  }

  if ( lsst->sol.cl[lsst->sol.nbcl].att == 'v' ) {
    for (i=0; i<lsst->info.dim; i++) {
      lsst->sol.cl[lsst->sol.nbcl].u[i] = u[i];
    }
  }
  else if ( lsst->sol.cl[lsst->sol.nbcl].att == 'n' ) {
    lsst->sol.cl[lsst->sol.nbcl].u[0] = u[0];
  }

	return(1);
}


/* specify gravity value */
void LS_setGra(LSst *lsst, double *gr) {
  int   i;

	lsst->info.load |= (1 << 0);
  for (i=0; i<lsst->info.dim; i++)
    lsst->info.gr[i] = gr[i];
}


/* specify elasticity Lame coefficients */
int LS_setLame(LSst *lsst,int ref,double lambda,double mu) {
  if ( lsst->sol.nmat == LS_MAT-1 )  return(0);

  lsst->sol.nmat++;
	lsst->sol.mat[lsst->sol.nmat].ref    = ref;
	lsst->sol.mat[lsst->sol.nmat].lambda = lambda;
	lsst->sol.mat[lsst->sol.nmat].mu     = mu;

  return(1);
}


/* construct mesh */
int LS_mesh(LSst *lsst,int np,int na,int nt,int ne) {
	int   dof;

	if ( !lsst )  return(0);

	lsst->info.np = np;
	lsst->info.na = na;
	lsst->info.nt = nt;
	lsst->info.ne = ne;

  /* bound on number of nodes */
  dof = lsst->info.typ == P2 ? 10 : 1;
  lsst->mesh.point = (pPoint)calloc(dof*lsst->info.np+1,sizeof(Point));
	assert(lsst->mesh.point);

  if ( lsst->info.na ) {
    lsst->mesh.edge  = (pEdge)calloc(lsst->info.na+1,sizeof(Edge));
    assert(lsst->mesh.edge);
  }
  if ( lsst->info.nt ) {
    lsst->mesh.tria  = (pTria)calloc(lsst->info.nt+1,sizeof(Tria));
    assert(lsst->mesh.tria);
  }
  if ( lsst->info.ne ) {
    lsst->mesh.tetra  = (pTetra)calloc(lsst->info.ne+1,sizeof(Tetra));
    assert(lsst->mesh.tetra);
  }

  return(1);
}

/* insert mesh elements into structure */
int LS_addVer(LSst *lsst,int idx,double *c,int ref) {
  pPoint   ppt;
	int      i;

	assert(idx > 0 && idx <= lsst->info.np);
	ppt = &lsst->mesh.point[idx];
	for (i=0; i<lsst->info.dim; i++)
	  ppt->c[i] = c[i];
	ppt->ref = ref;

	return(1);
}

int LS_addEdg(LSst *lsst,int idx,int *v,int ref) {
	pEdge   pe;
	
	assert(idx > 0 && idx <= lsst->info.na);
	pe = &lsst->mesh.edge[idx];
	memcpy(&pe->v[0],&v[0],2*sizeof(int));
  pe->ref = ref;

	return(1);
}

int LS_addTri(LSst *lsst,int idx,int *v,int ref) {
	pTria   pt;

	assert(idx > 0 && idx <= lsst->info.nt);
	pt = &lsst->mesh.tria[idx];
	memcpy(&pt->v[0],&v[0],3*sizeof(int));
  pt->ref = ref;

	return(1);
}

/* return mesh header */
void LS_headMesh(LSst *lsst,int *np,int *na,int *nt,int *ne) {
	*np = lsst->info.np;
	*na = lsst->info.na;
	*nt = lsst->info.nt;
	*ne = lsst->info.ne;
}

int LS_addTet(LSst *lsst,int idx,int *v,int ref) {
	pTria   pt;
	
	assert(idx > 0 && idx <= lsst->info.ne);
	pt = &lsst->mesh.tria[idx];
	memcpy(&pt->v[0],&v[0],4*sizeof(int));
  pt->ref = ref;

	return(1);
}


/* initialize solution vector or Dirichlet conditions 
   return: 1 if completion
           0 if no vertex array allocated
          -1 if previous data stored in struct. */
int LS_iniSol(LSst *lsst,double *u) {
  if ( !lsst->info.np )  return(0);

  /* no data already allocated */
  if ( !lsst->sol.u ) {
    lsst->sol.u  = (double*)u;
		return(1);
  }
	/* resolve potential conflict */
	else {
		free(lsst->sol.u);
    lsst->sol.u  = (double*)u;
		return(-1);
	}
}


/* return pointer to solution (Warning: starts at address 0) */
double *LS_getSol(LSst *lsst) {
	return(lsst->sol.u);
}


int LS_elastic(LSst *lsst) {
  int   ier;

  if ( lsst->info.dim == 2)
		ier = elasti1_2d(lsst);
	else
		ier = elasti1_3d(lsst);

	return(ier);	
}


