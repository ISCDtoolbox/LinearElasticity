#include "elastic.h"


/* compactify mesh structure */
int pack_3d(LSst *lsst) {
  pTetra    pei,pef;
  pTria     pti,ptf;
  double    l,m;
  int      *perm,i,k,nf;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pei = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pei->ref,&l,&m) ) {
      nf++;
      for (i=0; i<4; i++)  lsst->mesh.point[pei->v[i]].old = 1;
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

  /* compress and realloc vertices */
  lsst->info.npi = lsst->info.np;
  nf = 0;
  for (k=1; k<=lsst->info.np; k++) {
    if ( lsst->mesh.point[k].old > 0 ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
      lsst->mesh.point[nf].old = k;
      perm[k] = nf;
    }
  }
  for (k=nf+1; k<=lsst->info.np; k++)  lsst->mesh.point[k].old = 0;
  lsst->info.np = nf;

  /* compress and renum tetrahedra */
  lsst->info.nei = lsst->info.ne;
  nf = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pei = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pei->ref,&l,&m) ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.tetra[nf],&lsst->mesh.tetra[k],sizeof(Tetra));
      pef = &lsst->mesh.tetra[nf];
      for (i=0; i<4; i++)  pef->v[i] = perm[pef->v[i]];
    }
  }
  lsst->info.ne = nf;

  /* renum triangles */
  lsst->info.nti = lsst->info.nt;
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pti = &lsst->mesh.tria[k];
    if ( perm[pti->v[0]] && perm[pti->v[1]] && perm[pti->v[2]] ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
      ptf = &lsst->mesh.tria[nf];
      ptf->v[0] = perm[ptf->v[0]];
      ptf->v[1] = perm[ptf->v[1]];
      ptf->v[2] = perm[ptf->v[2]];
    }
  }
  lsst->info.nt = nf;

  /* compress solution (data) */
  if ( lsst->sol.u ) {
		for (k=1; k<=lsst->info.np; k++) {
			lsst->sol.u[3*(k-1)+0] = lsst->sol.u[3*(perm[k]-1)+0];
			lsst->sol.u[3*(k-1)+1] = lsst->sol.u[3*(perm[k]-1)+1];
			lsst->sol.u[3*(k-1)+2] = lsst->sol.u[3*(perm[k]-1)+2];
		}
  }
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    if ( lsst->info.ne )  fprintf(stdout,", %d tetrahedra",lsst->info.ne);
    fprintf(stdout,"\n");
  }
  free(perm);

  return(1);
}


int pack_2d(LSst *lsst) {
  pTria     pti,ptf;
  pEdge     pai,paf;
  double    l,m;
  int      *perm,i,k,nf;

  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pti = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pti->ref,&l,&m) ) {
      nf++;
      for (i=0; i<3; i++)  lsst->mesh.point[pti->v[i]].old = 1;
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
  lsst->info.npi = lsst->info.np;
  nf = 0;
  for (k=1; k<=lsst->info.np; k++) {
    if ( lsst->mesh.point[k].old > 0 ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
      lsst->mesh.point[nf].old = k;
      perm[k] = nf;
    }
  }
  for (k=nf+1; k<=lsst->info.np; k++)  lsst->mesh.point[k].old = 0;
  lsst->info.np = nf;

  /* compress and renum triangles */
  lsst->info.nti = lsst->info.nt;
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pti = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pti->ref,&l,&m) ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
      ptf = &lsst->mesh.tria[nf];
      for (i=0; i<3; i++)  ptf->v[i] = perm[ptf->v[i]];
    }
  }
  lsst->info.nt = nf;

  /* renum edges */
  nf = 0;
  for (k=1; k<=lsst->info.na; k++) {
    pai = &lsst->mesh.edge[k];
    if ( perm[pai->v[0]] && perm[pai->v[1]] ) {
      nf++;
      if ( nf < k )
        memcpy(&lsst->mesh.edge[nf],&lsst->mesh.edge[k],sizeof(Edge));
      paf = &lsst->mesh.edge[nf];
      paf->v[0] = perm[paf->v[0]];
      paf->v[1] = perm[paf->v[1]];
    }
  }
  lsst->info.na = nf;

  /* compress solution (data) */
  if ( lsst->sol.u ) {
		for (k=1; k<=lsst->info.np; k++) {
			lsst->sol.u[2*(k-1)+0] = lsst->sol.u[2*(perm[k]-1)+0];
			lsst->sol.u[2*(k-1)+1] = lsst->sol.u[2*(perm[k]-1)+1];
		}
  }

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    fprintf(stdout,"\n");
  }
  free(perm);

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
	unew = (double*)calloc(lsst->info.dim * (lsst->info.npi+lsst->info.np2),sizeof(double));
	assert(unew);

  for (k=1; k<=lsst->info.np; k++) {
    ppt = &lsst->mesh.point[k];
    if ( ppt->old ) {
			for (i=0; i<lsst->info.dim; i++)
			  unew[lsst->info.dim * (ppt->old-1)+i] = lsst->sol.u[lsst->info.dim * (k-1)+i];
    }
  }
	free(lsst->sol.u);
	lsst->sol.u   = unew;
  lsst->info.np = lsst->info.npi;
  lsst->info.na = 0;

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",lsst->info.np);
  }
  return(1);
}

