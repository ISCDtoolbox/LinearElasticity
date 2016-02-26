#include "elastic.h"
#include "ls_calls.h"
#include "libmesh5.h"


/* read mesh */
int loadMesh(LSst *lsst) {
  pPoint       ppt;
	pEdge        pa;
  pTria        pt1;
  pTetra       ptt;
  float        fp1,fp2,fp3;
  int          k,dof,inm;
  char        *ptr,data[256];

  strcpy(data,lsst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&lsst->info.ver,&lsst->info.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&lsst->info.ver,&lsst->info.dim)) ) {
        fprintf(stderr," # %s: file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&lsst->info.ver,&lsst->info.dim)) ) {
    fprintf(stderr," # %s: file not found.\n",data);
    return(0);
  }

  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  lsst->info.np = GmfStatKwd(inm,GmfVertices);
	lsst->info.na = GmfStatKwd(inm,GmfEdges);
  lsst->info.nt = GmfStatKwd(inm,GmfTriangles);
  lsst->info.ne = GmfStatKwd(inm,GmfTetrahedra);

  if ( !lsst->info.np ) {
    if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }
	lsst->info.npi = lsst->info.np;
	lsst->info.nai = lsst->info.na;
	lsst->info.nti = lsst->info.nt;
	lsst->info.nei = lsst->info.ne;

  /* memory alloc */
	dof = lsst->info.typ == P2 ? 10 : 1;  /* bound on number of nodes */
  lsst->mesh.point = (Point*)calloc(dof*lsst->info.np+1,sizeof(Point));
  assert(lsst->mesh.point);

  /* 2d mesh */
  GmfGotoKwd(inm,GmfVertices);
  if ( lsst->info.dim == 2 ) {
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      if ( lsst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->ref);
    }
  }
  else {
    for (k=1; k<=lsst->info.np; k++) {
      ppt = &lsst->mesh.point[k];
      if ( lsst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
  }
  if ( lsst->info.na > 0 ) {
    lsst->mesh.edge  = (pEdge)calloc(lsst->info.na+1,sizeof(Edge));
    assert(lsst->mesh.edge);
    /* read mesh edges */
    GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=lsst->info.na; k++) {
      pa = &lsst->mesh.edge[k];
      GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
    }
  }
  if ( lsst->info.nt > 0 ) {
    lsst->mesh.tria  = (Tria*)calloc(lsst->info.nt+1,sizeof(Tria));
    assert(lsst->mesh.tria);
    /* read mesh triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=lsst->info.nt; k++) {
      pt1 = &lsst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
  }
  if ( lsst->info.ne > 0 ) {
    lsst->mesh.tetra  = (Tetra*)calloc(lsst->info.ne+1,sizeof(Tetra));
    assert(lsst->mesh.tetra);
    /* read tetrahedra */
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=lsst->info.ne; k++) {
      ptt = &lsst->mesh.tetra[k];
      GmfGetLin(inm,GmfTetrahedra,&ptt->v[0],&ptt->v[1],&ptt->v[2],&ptt->v[3],&ptt->ref);
    }
  }
	GmfCloseMesh(inm);

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    if ( lsst->info.ne )  fprintf(stdout,", %d tetrahedra",lsst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load initial solution */
int loadSol(LSst *lsst) {
  float       buf[GmfMaxTyp];
  double      bufd[GmfMaxTyp];
  int         i,k,dim,ver,np,type,inm,typtab[GmfMaxTyp],offset;
  char       *ptr,data[128];

	if ( !lsst->sol.namein )  return(-1);
  strcpy(data,lsst->sol.namein);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != lsst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
  if ( !np || typtab[0] != GmfVec || np != lsst->info.np )  return(-1);

  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read mesh solutions */
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( ver == GmfFloat ) {
    for (k=0; k<lsst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,&buf);
    for (i=0; i<lsst->info.dim; i++)      
      lsst->sol.u[lsst->info.dim*k+i] = buf[i];
    }
  }
  else {
    for (k=0; k<lsst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
    for (i=0; i<lsst->info.dim; i++)      
        lsst->sol.u[lsst->info.dim*k+i] = bufd[i];
    }
  }
  GmfCloseMesh(inm);

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout," %d data vectors\n",lsst->info.np);
  }

  return(1);
}


int saveSol(LSst *lsst) {
  double    dbuf[GmfMaxTyp];
  float     fbuf[GmfMaxTyp];
  int       k,ia,i,outm,type,typtab[GmfMaxTyp];
  char     *ptr,data[128];

  strcpy(data,lsst->sol.nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,lsst->info.ver == 1 ? ".solb" : ".sol");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,lsst->info.ver,lsst->info.dim)) ) {
    fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);
  type = 1;
  typtab[0] = GmfVec;

  /* write sol */
  GmfSetKwd(outm,GmfSolAtVertices,lsst->info.np+lsst->info.np2,type,typtab);
  if ( lsst->info.ver == GmfFloat ) {
    for (k=0; k<lsst->info.np+lsst->info.np2; k++) {
      ia = lsst->info.dim * k;
      for (i=0; i<lsst->info.dim; i++)      
        fbuf[i] = lsst->sol.u[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
  }
  else {
    for (k=0; k<lsst->info.np+lsst->info.np2; k++) {
      ia = lsst->info.dim * k;
      for (i=0; i<lsst->info.dim; i++)      
        dbuf[i] = lsst->sol.u[ia+i];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }
  GmfCloseMesh(outm);

  if ( lsst->info.verb != '0' )  fprintf(stdout," %d data vectors\n",lsst->info.np+lsst->info.np2);

  return(1);
}


int saveMesh(LSst *lsst) {
  pPoint    ppt;
  pEdge     pe;
  pTria     ptt;
  pTetra    pt;
  double   *u;
  int       k,ia,i,outm;
  char     *ptr,data[128];

  strcpy(data,lsst->sol.nameout);
  ptr = strstr(data,".sol");
  if ( ptr )  *ptr = '\0';
  strcat(data,lsst->info.ver == 1 ? ".meshb" : ".mesh");

  lsst->info.ver = GmfDouble;
  if ( !(outm = GmfOpenMesh(data,GmfWrite,lsst->info.ver,lsst->info.dim)) ) {
    fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* write sol */
  GmfSetKwd(outm,GmfVertices,lsst->info.np);
  for (k=1; k<lsst->info.np; k++) {
    ppt = &lsst->mesh.point[k];
    u   = &lsst->sol.u[lsst->info.dim*(k-1)];
    if ( lsst->info.dim == 2 )
      GmfSetLin(outm,GmfVertices,ppt->c[0]+u[0],ppt->c[1]+u[1],ppt->ref);
    else
      GmfSetLin(outm,GmfVertices,ppt->c[0]+u[0],ppt->c[1]+u[1],ppt->c[2]+u[2],ppt->ref);
  }
  if ( lsst->info.nt > 0 ) {
    GmfSetKwd(outm,GmfTriangles,lsst->info.nt);
    for (k=1; k<=lsst->info.nt; k++) {
      ptt = &lsst->mesh.tria[k];
      GmfSetLin(outm,GmfTriangles,ptt->v[0],ptt->v[1],ptt->v[2],ptt->ref);
    }
  }
  if ( lsst->info.ne > 0 ) {
    GmfSetKwd(outm,GmfTetrahedra,lsst->info.ne);
    for (k=1; k<=lsst->info.ne; k++) {
      pt = &lsst->mesh.tetra[k];
      GmfSetLin(outm,GmfTetrahedra,pt->v[0],pt->v[1],pt->v[2],pt->v[3],pt->ref);
    }
  }
  if ( lsst->info.na > 0 ) {
    GmfSetKwd(outm,GmfEdges,lsst->info.na);
    for (k=1; k<=lsst->info.na; k++) {
      pe = &lsst->mesh.edge[k];
      GmfSetLin(outm,GmfEdges,pe->v[0],pe->v[1],pe->ref);
    }
  }
  GmfCloseMesh(outm);

  if ( lsst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",lsst->info.np);
    if ( lsst->info.na )  fprintf(stdout,", %d edges",lsst->info.na);
    if ( lsst->info.nt )  fprintf(stdout,", %d triangles",lsst->info.nt);
    if ( lsst->info.ne )  fprintf(stdout,", %d tetrahedra",lsst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1); 
} 