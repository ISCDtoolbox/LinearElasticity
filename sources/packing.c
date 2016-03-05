#include "elastic.h"


/* compactify mesh structure */
int pack_3d(LSst *lsst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    l,m,w[3];
  int       i,k,dof,nf,id;

  /* check if compression needed */
  nf  = 0;
  for (k=1; k<=lsst->info.ne; k++) {
    pe = &lsst->mesh.tetra[k];
    if ( getMat(&lsst->sol,pe->ref,&l,&m) ) {
      nf++;
      for (i=0; i<4; i++)  lsst->mesh.point[pe->v[i]].new = pe->v[i];
    }
  }
  if ( nf == lsst->info.nei )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;

  /* compress and renum vertices */
  nf = lsst->info.npi;
  k  = 1;
  while ( k <= nf ) {
    if ( lsst->mesh.point[k].new == 0 ) {
      while ( (lsst->mesh.point[nf].new == 0) && (k <= nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        /* swap solution too */
        if ( lsst->sol.u ) {
          memcpy(&w,&lsst->sol.u[3*(nf-1)],3*sizeof(double));
          memcpy(&lsst->sol.u[3*(nf-1)],&lsst->sol.u[3*(k-1)],3*sizeof(double));
          memcpy(&lsst->sol.u[3*(k-1)],&w,3*sizeof(double));
        }
        lsst->mesh.point[k].new  = nf;
        lsst->mesh.point[nf].new = k;
        nf--;
      }
    }
    k++;
  }
  lsst->info.np = nf;

  /* compress and renum tetrahedra */
  dof = lsst->info.typ == P1 ? 4 : 10;
  for (k=1; k<=lsst->info.nei; k++) {
    pe = &lsst->mesh.tetra[k];
    for (i=0; i<dof; i++)  pe->v[i] = lsst->mesh.point[pe->v[i]].new;
  }
  nf = lsst->info.nei;
  k  = 1;
  while ( k <= nf ) {
    pe = &lsst->mesh.tetra[k];
    if ( !getMat(&lsst->sol,pe->ref,&l,&m) ) {
      do {
        pe = &lsst->mesh.tetra[nf];
        if ( getMat(&lsst->sol,pe->ref,&l,&m) )  break;
        nf --;
      }
      while ( k <= nf );
      /* put nf into k */
      if ( k < nf ) {
        memcpy(&lsst->mesh.tetra[0],&lsst->mesh.tetra[nf],sizeof(Tetra));
        memcpy(&lsst->mesh.tetra[nf],&lsst->mesh.tetra[k],sizeof(Tetra));
        memcpy(&lsst->mesh.tetra[k],&lsst->mesh.tetra[0],sizeof(Tetra));
        nf--;
      }
    }
    k++;
  }
  lsst->info.ne = nf;

  /* renum triangles */
  dof = lsst->info.typ == P1 ? 3 : 6;
  for (k=1; k<=lsst->info.nti; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<dof; i++)  pt->v[i] = lsst->mesh.point[pt->v[i]].new;
  }
  nf = lsst->info.nti;
  k  = 1;
  while ( k <= nf ) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] > lsst->info.np )  break;
    if ( i < 3 ) {
      do {
        pt = &lsst->mesh.tria[nf];
        for (i=0; i<3; i++)
          if ( pt->v[i] > lsst->info.np )  break;
        if ( i == 3 )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&lsst->mesh.tria[0],&lsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
        memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[0],sizeof(Tria));
        nf--;
      }
    }
    k++;
  }
  lsst->info.nt = nf;

  /* renum edges */
  for (k=1; k<=lsst->info.nai; k++) {
    pa = &lsst->mesh.edge[k];
    for (i=0; i<3; i++)  pa->v[i] = lsst->mesh.point[pa->v[i]].new;
  }
  nf = lsst->info.nai;
  k  = 1;
  while ( k <= nf ) {
    pa = &lsst->mesh.edge[k];
    if ( (pa->v[0] == 0) || (pa->v[0] > lsst->info.np) || \
         (pa->v[1] == 0) || (pa->v[1] > lsst->info.np) ) {
      do {
        pa = &lsst->mesh.edge[nf];
        if ( (pa->v[0] > 0) && (pa->v[0] <= lsst->info.np) && \
             (pa->v[1] > 0) && (pa->v[1] <= lsst->info.np) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&lsst->mesh.edge[k],&lsst->mesh.edge[nf],sizeof(Edge));
        nf--;
      }
    }
    k++;
  }
  lsst->info.na = nf;

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
  int       i,k,dof,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=lsst->info.nti; k++) {
    pt = &lsst->mesh.tria[k];
    if ( getMat(&lsst->sol,pt->ref,&l,&m) ) {
      nf++;
      for (i=0; i<3; i++)  lsst->mesh.point[pt->v[i]].new = pt->v[i];
    }
  }
  if ( nf == lsst->info.nti )  return(-1);

  /* store permutations */
  if ( lsst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  lsst->info.zip = 1;

  /* compress and realloc vertices+solution/data */
  nf = lsst->info.npi;
  k  = 1;
  while ( k <= nf ) {
    if ( lsst->mesh.point[k].new == 0 ) {
      while ( (lsst->mesh.point[nf].new == 0) && (k <= nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&lsst->mesh.point[0],&lsst->mesh.point[nf],sizeof(Point));
        memcpy(&lsst->mesh.point[nf],&lsst->mesh.point[k],sizeof(Point));
        memcpy(&lsst->mesh.point[k],&lsst->mesh.point[0],sizeof(Point));
        /* swap solution too */
        if ( lsst->sol.u ) {
          memcpy(&w,&lsst->sol.u[2*(nf-1)],2*sizeof(double));
          memcpy(&lsst->sol.u[2*(nf-1)],&lsst->sol.u[2*(k-1)],2*sizeof(double));
          memcpy(&lsst->sol.u[2*(k-1)],&w,2*sizeof(double));
        }
        lsst->mesh.point[k].new  = nf;
        lsst->mesh.point[nf].new = k;
        nf--;
      }
    }
    k++;
  }
  lsst->info.np = nf;

  /* compress and renum triangles */
  dof = lsst->info.typ == P1 ? 3 : 6;
  for (k=1; k<=lsst->info.nti; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<dof; i++)  pt->v[i] = lsst->mesh.point[pt->v[i]].new;
  }
  nf = lsst->info.nti;
  k  = 1;
  while ( k <= nf ) {
    pt = &lsst->mesh.tria[k];
    if ( !getMat(&lsst->sol,pt->ref,&l,&m) ) {
      do {
        pt = &lsst->mesh.tria[nf];
        if ( getMat(&lsst->sol,pt->ref,&l,&m) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf into k */
      if ( k < nf ) {
        memcpy(&lsst->mesh.tria[0],&lsst->mesh.tria[nf],sizeof(Tria));
        memcpy(&lsst->mesh.tria[nf],&lsst->mesh.tria[k],sizeof(Tria));
        memcpy(&lsst->mesh.tria[k],&lsst->mesh.tria[0],sizeof(Tria));
        nf--;
      }
    }
    k++;
  }
  lsst->info.nt = nf;

  /* compress and renum edges */
  for (k=1; k<=lsst->info.nai; k++) {
    pa = &lsst->mesh.edge[k];
    for (i=0; i<3; i++)  pa->v[i] = lsst->mesh.point[pa->v[i]].new;
  }
  nf = lsst->info.nai;
  k  = 1;
  while ( k <= nf ) {
    pa = &lsst->mesh.edge[k];
    if ( (pa->v[0] == 0) || (pa->v[0] > lsst->info.np) || \
         (pa->v[1] == 0) ||(pa->v[1] > lsst->info.np) ) {
      do {
        pa = &lsst->mesh.edge[nf];
        if ( (pa->v[0] > 0) && (pa->v[0] <= lsst->info.np) && \
             (pa->v[1] > 0) && (pa->v[1] <= lsst->info.np) )  break;
        nf--;
      }
      while ( k <= nf );
      /* put nf in k */
      if ( k < nf ) {
        memcpy(&lsst->mesh.edge[k],&lsst->mesh.edge[nf],sizeof(Edge));
        nf--;
      }
    }
    k++;
  }
  lsst->info.na = nf;

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

