#include "elastic.h"


/* compactify mesh structure */
int pack_3d(LSst *lsst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    l,m,w[3];
  int       i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pe = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pe->ref,&l,&m) ) {
      nf++;
      for (i=0; i<4; i++)  lsst->mesh.point[pe->v[i]].new = pe->v[i];
    }
  }
  if ( nf == lsst->info.ne )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;

  /* compress and renum vertices */
  nf = lsst->info.ne;
  k  = 1;
  while ( k < nf ) {
    if ( lsst->mesh.point[k].new == 0 ) {
      while ( (lsst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        lsst->mesh.point[k].new  = nf;
        lsst->mesh.point[nf].new = k;
        
        if ( lsst->sol.u ) {
          memcpy(&w,&lsst->sol.u[3*(nf-1)],3*sizeof(double));
          memcpy(&lsst->sol.u[3*(nf-1)],&lsst->sol.u[3*(k-1)],3*sizeof(double));
          memcpy(&lsst->sol.u[3*(k-1)],&w,3*sizeof(double));
        }
      }
      nf++;
    }
    k++;
  }
  lsst->info.np = k;

  /* compress and renum tetrahedra */
  nf = lsst->info.ne;
  k  = 1;
  while ( k <= nf ) {
    pe = &lsst->mesh.tetra[k];
    if ( !getMat(&lsst->sol,pe->ref,&l,&m) ) {
      while ( !getMat(&lsst->sol,lsst->mesh.tetra[nf].ref,&l,&m) && (k < nf) )  nf --;
      /* put nf into k */
      memcpy(&lsst->mesh.tetra[k],&lsst->mesh.tetra[nf],sizeof(Tetra));
      nf--;
    }
    for (i=0; i<4; i++)  pe->v[i] = lsst->mesh.point[pe->v[i]].new;
    k++;
  }
  lsst->info.ne = k-1;

  /* renum triangles */
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    pt->v[0] = lsst->mesh.point[pt->v[0]].new;
    pt->v[1] = lsst->mesh.point[pt->v[1]].new;
    pt->v[2] = lsst->mesh.point[pt->v[2]].new;
  }
  
  /* renum edges */
  for (k=1; k<=lsst->info.na; k++) {
    pa = &lsst->mesh.edge[k];
    pa->v[0] = lsst->mesh.point[pa->v[0]].new;
    pa->v[1] = lsst->mesh.point[pa->v[1]].new;
  }

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    if ( lsst->info.ne )  fprintf(stdout,", %d tetrahedra",lsst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* mesh renumbering and packing */
int pack_2d(LSst *lsst) {
  pTria     pt;
  pEdge     pa;
  double    l,m,w[2];
  int       i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pt->ref,&l,&m) ) {
      nf++;
      for (i=0; i<3; i++)  lsst->mesh.point[pt->v[i]].new = pt->v[i];
    }
  }
  if ( nf == lsst->info.nt )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;

  /* compress and realloc vertices+solution/data */
  nf = lsst->info.np;
  k  = 1;
  while ( k < nf ) {
    if ( lsst->mesh.point[k].new == 0 ) {
      while ( (lsst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        lsst->mesh.point[k].new  = nf;
        lsst->mesh.point[nf].new = k;
        
        if ( lsst->sol.u ) {
          memcpy(&w,&lsst->sol.u[2*(nf-1)],2*sizeof(double));
          memcpy(&lsst->sol.u[2*(nf-1)],&lsst->sol.u[2*(k-1)],2*sizeof(double));
          memcpy(&lsst->sol.u[2*(k-1)],&w,2*sizeof(double));
        }
      }
      nf--;
    }
    k++;
  }
  lsst->info.np = k;

  /* compress and renum triangles */
  nf = lsst->info.nt;
  k  = 1;
  while ( k <= nf ) {
    pt = &lsst->mesh.tria[k];
    if ( !getMat(&lsst->sol,pt->ref,&l,&m) ) {
      while ( !getMat(&lsst->sol,lsst->mesh.tria[nf].ref,&l,&m) && (k < nf) )  nf --;
      /* put nf into k */
      memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[nf],sizeof(Tria));
      nf--;
    }
    for (i=0; i<3; i++)  pt->v[i] = lsst->mesh.point[pt->v[i]].new;
    k++;
  }
  lsst->info.nt = k-1;

  /* simply renum edges */
  for (k=1; k<=lsst->info.na; k++) {
    pa = &lsst->mesh.edge[k];
    pa->v[0] = lsst->mesh.point[pa->v[0]].new;
    pa->v[1] = lsst->mesh.point[pa->v[1]].new;
  }

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
	double  w[3];
  int     k,dim;
	char    i;

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Uncompressing data: ");
    fflush(stdout);
  }

  dim = lsst->info.dim;
  for (k=1; k<=lsst->info.np; k++) {
    ppt = &lsst->mesh.point[k];
    if ( ppt->new != k ) {
      memcpy(&w,&lsst->sol.u[dim*(k-1)+0],dim*sizeof(double));
      memcpy(&lsst->sol.u[dim*(k-1)+0],&lsst->sol.u[dim*(ppt->new-1)+0],dim*sizeof(double));
      memcpy(&lsst->sol.u[dim*(ppt->new-1)+0],&w,dim*sizeof(double));
    }
  }
  lsst->info.np = lsst->info.npi;
  lsst->info.na = lsst->info.nai;
  lsst->info.nt = lsst->info.nti;
  lsst->info.ne = lsst->info.nei;

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",lsst->info.np);
  }

  return(1);
}

