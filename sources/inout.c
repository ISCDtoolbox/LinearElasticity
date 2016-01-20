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
        fprintf(stderr," # %s file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&lsst->info.ver,&lsst->info.dim)) ) {
    fprintf(stderr," # %s file not found.\n",data);
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
	lsst->info.nti = lsst->info.nt;
	lsst->info.nei = lsst->info.ne;

  /* memory alloc */
	dof = lsst->info.typ == P2 ? 10 : 1;  /* bound on number of nodes */
  lsst->mesh.point = (Point*)calloc(dof*lsst->info.np+1,sizeof(Point));
  assert(lsst->mesh.point);
  if ( lsst->info.nt > 0 ) {
    lsst->mesh.tria  = (Tria*)calloc(lsst->info.nt+1,sizeof(Tria));
    assert(lsst->mesh.tria);
  }
  if ( lsst->info.ne ) {
    lsst->mesh.tetra  = (Tetra*)calloc(lsst->info.ne+1,sizeof(Tetra));
    assert(lsst->mesh.tetra);
  }

  /* 2d mesh */
  if ( lsst->info.dim == 2 ) {
    GmfGotoKwd(inm,GmfVertices);
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
    /* read mesh edges */
    if ( lsst->info.na ) {
      lsst->mesh.edge  = (pEdge)calloc(lsst->info.na+1,sizeof(Edge));
      assert(lsst->mesh.edge);
    }
    GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=lsst->info.na; k++) {
      pa = &lsst->mesh.edge[k];
      GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
    }
    /* read mesh triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=lsst->info.nt; k++) {
      pt1 = &lsst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
  }
	/* 3d mesh */
  else {
    GmfGotoKwd(inm,GmfVertices);
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
    /* read triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=lsst->info.nt; k++) {
      pt1 = &lsst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
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
    strcat(data,".solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".sol");
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,lsst->info.ver,lsst->info.dim)) ) {
    fprintf(stderr," # nable to open %s\n",data);
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
