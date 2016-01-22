#include "elastic.h"


/* compactify mesh structure */
int pack_3d(LSst *lsst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    l,m;
  int      *perm,i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pe = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pe->ref,&l,&m) ) {
      nf++;
      for (i=0; i<4; i++)  lsst->mesh.point[pe->v[i]].old = pe->v[i];
    }
  }
  if ( nf == lsst->info.ne )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;
  perm = (int*)calloc(lsst->info.np+1,sizeof(int));
  assert(perm);

  /* compress and renum vertices */
  nf = 0;
  for (k=1; k<=lsst->info.np; k++) {
    perm[k] = k;
    id = lsst->mesh.point[k].old;
    lsst->mesh.point[k].old = k;
    if ( id > 0 ) {
      nf++;
      /* swap k and nf */
      if ( nf < k ) {
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        perm[k] = nf;
      }
    }
  }
  lsst->info.np = nf;

  /* renum edges */
  nf = 0;
  for (k=1; k<=lsst->info.na; k++) {
    pa = &lsst->mesh.edge[k];
    if ( perm[pa->v[0]] <= lsst->info.np && perm[pa->v[1]] <= lsst->info.np ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.edge[nf],&lsst->mesh.edge[k],sizeof(Edge));
      pa = &lsst->mesh.edge[nf];
      pa->v[0] = perm[pa->v[0]];
      pa->v[1] = perm[pa->v[1]];
    }
  }
  lsst->info.na = nf;

  /* compress and renum triangles */
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pt->ref,&l,&m) ) {
      nf++;
      if ( nf < k ) {
        memcpy(&lsst->mesh.tria[0],&lsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
        memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[0],sizeof(Tria));
      }
    }
  }
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = perm[pt->v[i]];
  }
  lsst->info.nt = nf;

  /* compress and renum tetrahedra */
  nf = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pe = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pe->ref,&l,&m) ) {
      nf++;
      if ( nf < k ) {
        memcpy(&lsst->mesh.tria[0],&lsst->mesh.tria[nf],sizeof(Tetra));
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tetra));
        memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[0],sizeof(Tetra));
      }
    }
  }
  for (k=1; k<=lsst->info.ne; k++) {
    pe = &lsst->mesh.tetra[k];
    for (i=0; i<4; i++)  pe->v[i] = perm[pe->v[i]];
  }
  lsst->info.ne = nf;

  /* compress solution (data) */
  if ( lsst->sol.u ) {
		for (k=1; k<=lsst->info.npi; k++) {
			lsst->sol.u[2*(perm[k]-1)+0] = lsst->sol.u[2*(k-1)+0];
			lsst->sol.u[2*(perm[k]-1)+1] = lsst->sol.u[2*(k-1)+1];
			lsst->sol.u[2*(perm[k]-1)+2] = lsst->sol.u[2*(k-1)+2];
		}
  }
  free(perm);

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    if ( lsst->info.ne )  fprintf(stdout,", %d tetrahedra",lsst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


int pack_2d(LSst *lsst) {
  pTria     pt;
  pEdge     pa;
  double    l,m;
  int      *perm,i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pt->ref,&l,&m) ) {
      nf++;
      for (i=0; i<3; i++)  lsst->mesh.point[pt->v[i]].old = pt->v[i];
    }
  }
  if ( nf == lsst->info.nt )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;
  perm = (int*)calloc(lsst->info.np+1,sizeof(int));
  assert(perm);

  /* compress and realloc vertices */
  nf = 0;
  for (k=1; k<=lsst->info.np; k++) {
    perm[k] = k;
    id = lsst->mesh.point[k].old;
    lsst->mesh.point[k].old = k;
    if ( id > 0 ) {
      nf++;
      /* swap k and nf */
      if ( nf < k ) {
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        perm[k] = nf;
      }
    }
  }
  lsst->info.np = nf;

  /* compress and renum triangles */
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pt->ref,&l,&m) ) {
      nf++;
      if ( nf < k ) {
        memcpy(&lsst->mesh.tria[0],&lsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
        memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[0],sizeof(Tria));
      }
    }
  }
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = perm[pt->v[i]];
  }
  lsst->info.nt = nf;

  /* renum edges */
  nf = 0;
  for (k=1; k<=lsst->info.na; k++) {
    pa = &lsst->mesh.edge[k];
    if ( perm[pa->v[0]] <= lsst->info.np && perm[pa->v[1]] <= lsst->info.np ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.edge[nf],&lsst->mesh.edge[k],sizeof(Edge));
      pa = &lsst->mesh.edge[nf];
      pa->v[0] = perm[pa->v[0]];
      pa->v[1] = perm[pa->v[1]];
    }
  }
  lsst->info.na = nf;

  /* compress solution (data) */
  if ( lsst->sol.u ) {
		for (k=1; k<=lsst->info.npi; k++) {
			lsst->sol.u[2*(perm[k]-1)+0] = lsst->sol.u[2*(k-1)+0];
			lsst->sol.u[2*(perm[k]-1)+1] = lsst->sol.u[2*(k-1)+1];
		}
  }
  free(perm);

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* restore solution at initial vertices */
int unpack(LSst *lsst) {
  pPoint  ppt;
	double *unew;
  int     k;
	char    i;

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Uncompressing data: ");
    fflush(stdout);
  }
	unew = (double*)calloc(lsst->info.dim*lsst->info.npi,sizeof(double));
	assert(unew);

  for (k=1; k<=lsst->info.npi; k++) {
    ppt = &lsst->mesh.point[k];
		for (i=0; i<lsst->info.dim; i++)
			unew[lsst->info.dim*(ppt->old-1)+i] = lsst->sol.u[lsst->info.dim*(k-1)+i];
  }
	free(lsst->sol.u);
	lsst->sol.u   = unew;
  lsst->info.np = lsst->info.npi;
  lsst->info.na = lsst->info.nai;
  lsst->info.nt = lsst->info.nti;
  lsst->info.ne = lsst->info.nei;

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",lsst->info.np);
  }
  return(1);
}

