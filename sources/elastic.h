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
#define LS_REL   "Sep.2, 2015"
#define LS_CPY   "(C) Copyright 2010-2015, ICS-SU"
#define LS_STR   "------------------------------------------------------------"

#define LS_LAMBDA     10.0e5
#define LS_MU          8.2e5
#define LS_MAT        50
#define LS_CL         50
#define LS_RES        1.0e-6
#define LS_MAXIT      10000
#define LS_TGV        1.e+30

#define LS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define LS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )

#define EPS      1.e-6
#define EPSD     1.e-30


typedef struct {
  double   c[3];
  int      ref,old;
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
	int      np,np2,na,nt,ne,npi,nti,nei;
  double   gravity,gr[3];
  char     imprim,ddebug,load,typ,cg,zip,rhs,mfree;
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
  int      nit,iter,nbcl,nmat;
  double  *u,*bc,err;
  char    *namein,*nameout,cltyp;
  Cl      *cl;
  Mat     *mat;
} Sol;
typedef Sol * pSol;

/* hashing structure */
typedef struct {
  int   ia,ib,k,nxt;
} hedge;

typedef struct {
  int     siz,max,nxt;
  hedge  *item;
} Hash;

struct _LSst{
  Mesh    mesh;
	Sol     sol;
	Info    info;
	Hash    hash;
};

/* prototypes */
int  loadMesh(LSst *lsst);
int  loadSol(LSst *lsst);
int  saveSol(LSst *lsst);
int  pack_2d(LSst *lsst);
int  pack_3d(LSst *lsst);
int  unpack(LSst *lsst);
int  hashar_2d(LSst *lsst);
int  hashar_3d(LSst *lsst);
int  hashP2(Hash *hash,int a,int b);
pCl  getCl(pSol sol,int ref,int elt);
int  getMat(pSol sol,int ref,double *lambda,double *mu);

/* function pointers */
int   elasti1_2d(LSst *lsst);
int   elasti1_3d(LSst *lsst);



#endif