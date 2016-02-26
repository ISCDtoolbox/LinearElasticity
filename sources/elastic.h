#ifndef _ELASTIC_H
#define _ELASTIC_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "ls_calls.h"

#define LS_VER   "5.0c"
#define LS_REL   "Jan.19, 2016"
#define LS_CPY   "(C) Copyright 2006- , ICS-SU"

#define LS_LAMBDA     10.0e5
#define LS_MU          8.2e5
#define LS_MAT        50
#define LS_CL         50
#define LS_RES        1.0e-6
#define LS_MAXIT      10000
#define LS_TGV        1.e+30
#define LS_EPSD       1.e-200

#define LS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define LS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )


/* data structures */
typedef struct {
  double   c[3];
  int      ref,new;
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[3],ref;
} Edge;
typedef Edge * pEdge;

typedef struct {
  int     v[6],ref;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[10],ref;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
	int      dim,ver;
	int      np,np2,na,nt,ne,npi,nai,nti,nei;
  char     verb,typ,zip,mfree,xport;
  mytime   ctim[TIMEMAX];
} Info;

typedef struct {
  char    *name;
  pPoint   point;
	pEdge    edge;
  pTria    tria;
  pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  double   u[3];
  int      ref;
  char     typ,elt,att;
} Cl;
typedef Cl * pCl;

typedef struct {
  double  lambda,mu;
  int     ref;
} Mat;
typedef Mat * pMat;

typedef struct {
  int      nit,nbcl,nmat;
  double  *u,*F,res,gr[3];
  char    *namein,*nameout,*namepar,cltyp,clelt;
  Cl      *cl;
  Mat     *mat;
} Sol;
typedef Sol * pSol;

struct _LSst {
  Mesh    mesh;
	Sol     sol;
	Info    info;
};

/* prototypes */
int  loadMesh(LSst *lsst);
int  loadSol(LSst *lsst);
int  saveSol(LSst *lsst);
int  saveMesh(LSst *lsst);
int  pack_2d(LSst *lsst);
int  pack_3d(LSst *lsst);
int  unpack(LSst *lsst);
int  hashar_2d(LSst *lsst);
int  hashar_3d(LSst *lsst);
pCl  getCl(pSol sol,int ref,int elt);
int  getMat(pSol sol,int ref,double *lambda,double *mu);
int  elasti1_2d(LSst *lsst);
int  elasti1_3d(LSst *lsst);


#endif
