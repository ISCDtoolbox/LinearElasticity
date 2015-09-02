#ifndef _SM_H
#define _SM_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "lplib3.h"

#define CS_VER   "2.0a"
#define CS_REL   "Feb 18, 2012"
#define CS_CPY   "Copyright (c) LJLL, 2006-2012"

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define CS_TGV2    1.e+60
#define CS_NUL     1.e-150

#define CS_MAXCPU  128
#define CS_MAXIT   99999
#define CS_DIAG    (1<<0)
#define CS_UT      (1<<1)
#define CS_LT      (1<<2)
#define CS_SYM     (1<<3)

/* compressed sparse row (CSR) format */
typedef struct {
  double    *val;
  int       *col,*row,nr,nc,nbe,nmax;
	char       typ;
	void      *hm;
} Csr;
typedef Csr * pCsr;

/* compressed sparse vector (CSV) format */
typedef struct {
	double    *val;
	int       *col,nr,nc,nbe;
	void      *hm;
} Csv;
typedef Csv * pCsv;


/* for multithreading */
typedef struct {
	Csr    *A;
	double *x,*y,*z,l,m,r[CS_MAXCPU];
} CsrArg;


#define CS_MAT(M) (M && (M->nbe >= 0))
#define CS_CSR(M) (M && (M->nbe == -1))

void csrInit(int nbcpu);
void csrStop();

pCsv csvNew(int nr,int nc,int nmax);
int  csvFree(pCsv M);
int  csvPack(pCsv M);
int  csvPut(pCsv M,int i,int j,double val);
int  csvSet(pCsv M,int i,int j,double val);

pCsr csrNew(int nr,int nc,int nmax,char typ);
void csrAlloc(pCsr M,int nr,int nc,int nmax,char typ);
int  csrFree(pCsr M);
int  csrPack(pCsr M);
int  csrPut(pCsr M,int i,int j,double val);
int  csrSet(pCsr M,int i,int j,double val);
int  csrGet(pCsr M,int i,int j,double *val);
pCsr csrTr(pCsr M);
pCsr csrMulAAt(pCsr A);
pCsr csrMul(pCsr A,pCsr B);
pCsr csrAdd(pCsr A,pCsr B,double l,double m);
int  csrAx(pCsr A,double *x,double *y);
int  csrAtx(pCsr A,double *x,double *y);
int  csrAxpy(pCsr A,double *x,double *y,double *z,double l,double m); 
int  csrAtxpy(pCsr A,double *x,double *y,double *z,double l,double m); 
void csrlX(double *x,double *y,double l,int n);
void csrlXmY(double *x,double *y,double *z,double l,double m,int n);
double csrAxdotx(pCsr A,double *x,double *y);
double csrXY(double *x,double *y,int n);
double csrNorm(pCsr M);

int  csrPrint(pCsr M,int brief);
void csrPrLine(pCsr A,int i);
void csrPrVal(pCsr A,int i,int j);
pCsr csrLoad(char *name);
void csrSave(pCsr A,char *name);


/* solve */
int  csrSSOR(pCsr A,pCsr L,double *x,double *b);
int  csrGradient(pCsr A,double *x,double *b,double *err,int *nit);
int  csrConjGrad(pCsr A,double *x,double *b,double *err,int *nit);
int  csrConjGradGen(pCsr A,double *x,double *b,double *,char *,char ,double *err,int *nit);
int  csrPrecondGrad(pCsr A,double *x,double *b,double *er,int *ni,char tgv);
int  csrGMRES(pCsr A,double *x,double *b,double *er,int *ni,int krylov,int prec);


#endif
