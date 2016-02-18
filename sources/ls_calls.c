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
  lsst->sol.res = LS_RES;
  lsst->sol.nit = LS_MAXIT;
  lsst->sol.nbcl = 0;
  lsst->sol.nmat = 0;

  /* global parameters */
  lsst->info.dim    = dim;
  lsst->info.ver    = ver;
  lsst->info.verb   = '1';
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
  free(lsst->sol.u);
	free(lsst->sol.cl);
	free(lsst->sol.mat);

  chrono(OFF,&lsst->info.ctim[0]);
  if ( lsst->info.verb != '0' ) {
	  printim(lsst->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s sec.\n",stim);
  }

	return(1);
}

/* set params (facultative): verb= '-|0|+',  zip = 0|1 */
void LS_setPar(LSst *lsst,char verb,int zip) {
	lsst->info.verb = verb;
	lsst->info.zip  = zip;
}

/* handle boundary conditions:
  typ= Dirichlet, Load
  ref= integer
  att= char 'v', 'f', 'n'
  elt= enum LS_ver, LS_edg, LS_tri, LS_tet */
int LS_setBC(LSst *lsst,int typ,int ref,char att,int elt,double *u) {
  Cl    *pcl;
  int    i;

  pcl = &lsst->sol.cl[lsst->sol.nbcl];
  pcl->typ = typ;
  pcl->ref = ref;
  pcl->att = att;
  pcl->elt = elt;

  if ( pcl->typ == Dirichlet ) {
    if ( !strchr("fv",pcl->att) ) {
      fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
  }
  else if ( pcl->typ == Load ) {
    if ( !strchr("fnv",pcl->att) ) {
      if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
    if ( (pcl->elt == LS_ver) && (pcl->att == 'n') ) {
      if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: %c\n",pcl->att);
      return(0);
    }
  }

  if ( pcl->att == 'v' ) {
    for (i=0; i<lsst->info.dim; i++)  pcl->u[i] = u[i];
  }
  else if ( pcl->att == 'n' ) {
    pcl->u[0] = u[0];
  }

  if ( lsst->sol.nbcl == LS_CL-1 )  return(0);
  lsst->sol.nbcl++;

  return(1);
}


/* specify gravity value */
void LS_setGra(LSst *lsst, double *gr) {
  int   i;

	lsst->sol.cltyp |= Gravity;
  for (i=0; i<lsst->info.dim; i++)
    lsst->sol.gr[i] = gr[i];
}


/* specify elasticity Lame coefficients */
int LS_setLame(LSst *lsst,int ref,double lambda,double mu) {
  Mat    *pm;

  if ( lsst->sol.nmat == LS_MAT-1 )  return(0);
  
  pm = &lsst->sol.mat[lsst->sol.nmat];
  pm->ref    = ref;
  pm->lambda = lambda;
  pm->mu     = mu;
  
  lsst->sol.nmat++;

  return(1);
}

/* Construct solution */
int LS_newSol(LSst *lsst) {
  
  lsst->sol.u = (double*)calloc(lsst->info.dim * (lsst->info.np+lsst->info.np2),sizeof(double));
  assert(lsst->sol.u);
  return(1);
}

/* Add element u[dim] to solution at ip */
int LS_addSol(LSst *lsst,int ip,double *u) {
  memcpy(&lsst->sol.u[lsst->info.dim*(ip-1)],u,lsst->info.dim*sizeof(double));
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

/* store all mesh vertices + references */
int LS_allVer(LSst *lsst,int np,double *c,int *ref) {
  pPoint   ppt;
  int      i,k,dof;

	if ( !lsst )  return(0);

	lsst->info.np = np;
  dof = lsst->info.typ == P2 ? 10 : 1;
  lsst->mesh.point = (pPoint)calloc(dof*lsst->info.np+1,sizeof(Point));
	assert(lsst->mesh.point);

  for (k=1; k<=lsst->info.np; k++) {
    ppt = &lsst->mesh.point[k];
  	for (i=0; i<lsst->info.dim; i++)
      ppt->c[i] = c[lsst->info.dim*(k-1)+1];
  	ppt->ref = ref[k];
  }

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

/* store all mesh edges + references */
int LS_allEdg(LSst *lsst,int na,int *edg,int *ref) {
  pEdge    pe;
  int      k;

	if ( !lsst )  return(0);

  lsst->info.na = na;
  if ( lsst->info.na ) {
    lsst->mesh.edge  = (pEdge)calloc(lsst->info.na+1,sizeof(Edge));
    assert(lsst->mesh.edge);
  }

  for (k=1; k<=lsst->info.na; k++) {
    pe = &lsst->mesh.edge[k];
    pe->v[0] = edg[2*(k-1)+1];
  	pe->ref = ref[k];
  }

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

/* store all mesh edges + references */
int LS_allTri(LSst *lsst,int nt,int *tri,int *ref) {
  pTria    pt;
  int      i,k;

	if ( !lsst )  return(0);

  lsst->info.nt = nt;
  if ( lsst->info.nt ) {
    lsst->mesh.tria  = (pTria)calloc(lsst->info.nt+1,sizeof(Tria));
    assert(lsst->mesh.tria);
  }

  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<3; i++)
      pt->v[i] = tri[3*(k-1)+i+1];
    pt->ref = ref[k];
  }

  return(1);
}

int LS_addTet(LSst *lsst,int idx,int *v,int ref) {
	pTetra   pt;

	assert(idx > 0 && idx <= lsst->info.ne);
	pt = &lsst->mesh.tetra[idx];
	memcpy(&pt->v[0],&v[0],4*sizeof(int));
  pt->ref = ref;

	return(1);
}

/* store all mesh edges + references */
int LS_allTet(LSst *lsst,int ne,int *tet,int *ref) {
  pTetra   pt;
  int      i,k;

	if ( !lsst )  return(0);

  lsst->info.ne = ne;
  if ( lsst->info.ne ) {
    lsst->mesh.tetra  = (pTetra)calloc(lsst->info.ne+1,sizeof(Tetra));
    assert(lsst->mesh.tetra);
  }

  for (k=1; k<=lsst->info.ne; k++) {
    pt = &lsst->mesh.tetra[k];
    for (i=0; i<4; i++)
      pt->v[i] = tet[4*(k-1)+i+1];
    pt->ref = ref[k];
  }

  return(1);
}

/* return mesh header */
void LS_headMesh(LSst *lsst,int *np,int *na,int *nt,int *ne) {
  int k;
	*np = lsst->info.np;
	*na = lsst->info.na;
	*nt = lsst->info.nt;
	*ne = lsst->info.ne;
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

/* initialize right-hand side */
int LS_iniRHS(LSst *lsst,double *F) {
  /* no data already allocated */
  if ( !lsst->sol.F ) {
    lsst->sol.F  = (double*)F;
		return(1);
  }
	/* resolve potential conflict */
	else {
		free(lsst->sol.F);
    lsst->sol.F  = (double*)F;
		return(-1);
	}
}

/* return pointer to solution (Warning: starts at address 0) */
double *LS_getSol(LSst *lsst) {
	return(lsst->sol.u);
}


int LS_elastic(LSst *lsst) {
  int   i,ier;
  Cl    *pcl;

  for (i=0; i<lsst->sol.nbcl; i++) {
    pcl = &lsst->sol.cl[i];
		lsst->sol.cltyp |= pcl->typ;
    lsst->sol.clelt |= pcl->elt;
  }

  if ( lsst->info.dim == 2)
		ier = elasti1_2d(lsst);
	else
		ier = elasti1_3d(lsst);

	return(ier);	
}


