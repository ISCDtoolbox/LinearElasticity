#include "sparse.h"

static int CSR_libId  = 0;
static int CSR_libCpu = 1;

/* hashing for matrix storage */
typedef struct {
  double  val;
  int     j,nxt;
} hmat;

typedef struct {
  int    siz,max,nxt,nbe;
  hmat  *item;
} csHash;

/* define hash table for matrix */
static csHash *cshNew(int nr,int hmax) {
  csHash  *hm;
  int      k;

  hm = (csHash*)calloc(1,sizeof(csHash));
  assert(hm);

  /* adjust hash table params */
  hm->item = (hmat*)calloc(hmax+1,sizeof(hmat));
  assert(hm->item);
  hm->siz  = nr;
  hm->nxt  = nr;
  hm->max  = hmax + 1;
  hm->nbe  = 0;
  for (k=0; k<hm->siz; k++)  hm->item[k].j = -1;
  for (k=hm->siz; k<hm->max; k++)  hm->item[k].nxt = k+1;

  return(hm);
}

/* add val to a(i,j) via hash table */
static int cshPut(csHash *hm,int i,int j,double val) {
  hmat   *phm;
  int     key;

  key = i % hm->siz;
  phm = &hm->item[key];

  if ( phm->j == j ) {
    phm->val += val;
    return(1);
  }
  else if ( phm->j > -1 ) {
    while ( phm->nxt && phm->nxt < hm->max ) {
      phm = &hm->item[phm->nxt];
      if ( phm->j == j ) {
        phm->val += val;
        return(1);
      }
    }
    phm->nxt = hm->nxt;
    phm      = &hm->item[hm->nxt];
    phm->j   = j;
    phm->val = val;
    phm->nxt = 0;
    ++hm->nxt;
    hm->nbe++;
    /* check for overflow */
    if ( hm->nxt >= hm->max ) {
      /*fprintf(stdout,"    realloc de %d a %d\n",hm->max,(int)(1.4*hm->max));*/
      hm->max *= 1.4;
      hm->item = (hmat*)realloc(hm->item,hm->max*sizeof(hmat));
      assert(hm->item);
      for (j=hm->nxt; j<hm->max; j++)  hm->item[j].nxt = j+1;
    }
    return(1);
  }
  hm->nbe++;
  phm->j   = j;
  phm->val = val;
  phm->nxt = 0;
  return(1);
}

/* set val to a(i,j) via hash table */
static int cshSet(csHash *hm,int i,int j,double val) {
  hmat   *phm;
  int     key;

  key = i % hm->siz;
  phm = &hm->item[key];
  if ( phm->j == j ) {
    phm->val = val;
    return(1);
  }
  else if ( phm->j > -1 ) {
    while ( phm->nxt && phm->nxt < hm->max ) {
      phm = &hm->item[phm->nxt];
      if ( phm->j == j ) {
        phm->val = val;
        return(1);
      }
    }
    phm->nxt = hm->nxt;
    phm      = &hm->item[hm->nxt];
    phm->j   = j;
    phm->val = val;
    phm->nxt = 0;
    ++hm->nxt;
    hm->nbe++;
    /* check for overflow */
    if ( hm->nxt >= hm->max ) {
      /*printf("    realloc de %d a %d\n",hm->max,(int)(1.4*hm->max));*/
      hm->max *= 1.4;
      hm->item  = (hmat*)realloc(hm->item,hm->max*sizeof(hmat));
      assert(hm->item);
      for (j=hm->nxt; j<hm->max; j++)  hm->item[j].nxt = j+1;
    }
    return(1);
  }
  hm->nbe++;
  phm->j   = j;
  phm->val = val;
  phm->nxt = 0;
  return(1);
}

static pCsr csmNew(int nr,int nc,int nbe,char typ) {
  pCsr     C;

  C = (Csr*)calloc(1,sizeof(Csr));
  C->typ = typ;
  C->nr  = nr;
  C->nc  = nc;
  C->nbe = C->nmax = nbe;   /* upper bound */
  C->row = (int*)malloc((C->nr+1)*sizeof(int));
  C->col = (int*)malloc(C->nbe*sizeof(int));
  C->val = (double*)malloc(C->nbe*sizeof(double));
  assert(C->row);
  assert(C->col);
  assert(C->val);

  return(C);
}

/* set val to a(i,j) via matrix structure */
static int csmSet(pCsr M,int i,int j,double v) {
  int     k;

  if ( i >= M->nr || j >= M->nc )  return(0);
  for (k=M->row[i]; k<M->row[i+1]; k++) {
    if ( M->col[k] == j ) {
      M->val[k] = v;
      return(1);
    }
  }
  /* need realloc structure */
  M->nbe++;
  if ( M->nbe >= M->nmax ) {
    M->nmax *= 1.4;
    M->col = realloc(M->col,M->nmax*sizeof(int));
    assert(M->col);
    M->val = realloc(M->val,M->nmax*sizeof(double));
    assert(M->val);
  }
  for (k=i+1; k<M->nr; k++)  M->row[k] += 1;
  for (k=M->nbe-1; k>M->row[i]; k--)  M->col[k] = M->col[k-1];
  for (k=M->nbe-1; k>M->row[i]; k--)  M->val[k] = M->val[k-1];

  /* insert new cell */
  M->col[k] = j;
  M->val[k] = v;
  return(1);
}

/* add val to a(i,j) via matrix structure */
static int csmPut(pCsr M,int i,int j,double v) {
  int     k;

  if ( i >= M->nr || j >= M->nc )  return(0);
  for (k=M->row[i]; k<M->row[i+1]; k++) {
    if ( M->col[k] == j ) {
      M->val[k] += v;
      return(1);
    }
  }
  /* need realloc structure: but inefficient... */
  M->nbe++;
  if ( M->nbe >= M->nmax ) {
    M->nmax *= 1.4;
    M->col = realloc(M->col,M->nmax*sizeof(int));
    assert(M->col);
    M->val = realloc(M->val,M->nmax*sizeof(double));
    assert(M->val);
  }
  for (k=i+1; k<M->nr; k++)  M->row[k] += 1;
  for (k=M->nbe-1; k>M->row[i]; k--)  M->col[k] = M->col[k-1];
  for (k=M->nbe-1; k>M->row[i]; k--)  M->val[k] = M->val[k-1];

  /* insert new cell */
  M->col[k] = j;
  M->val[k] = v;
  return(1);
}

/* return a(i,j) if stored, 0 otherwise */
int csrGet(pCsr M,int i,int j,double *val) {
  int   k;

  *val = 0.0;
  if ( i >= M->nr || j >= M->nc )   return(0);
  for (k=M->row[i]; k<M->row[i+1]; k++) {
    if ( M->col[k] == j ) {
      *val = M->val[k];
      return(1);
    }
  }
  return(0);
}

/* get a(i,j) if stored */
static int cshGet(csHash *hm,int i,int j,double *val) {
  hmat   *phm;
  int     key;

  key = i % hm->siz;
  phm = &hm->item[key];
  *val = 0.0;

  if ( phm->j == j ) {
    *val = phm->val;
    return(1);
  }
  else {
    while ( phm->nxt && phm->nxt < hm->max ) {
      phm = &hm->item[phm->nxt];
      if ( phm->j == j ) {
        *val = phm->val;
        return(1);
      }
    }
    return(0);
  }
  return(0);
}


/*======= Part 1: CSV matrix =======*/

/* alloc memory to store matrix */
pCsv csvNew(int nr,int nc,int nmax) {
  pCsv   M;

  M = (pCsv)calloc(1,sizeof(Csr));  
  assert(M);
  M->nr   = nr;
  M->nc   = nc;
  M->nbe  = 0;
  M->hm   = cshNew(nr,nmax);

  return(M);
}

/* release memory */
int csvFree(pCsv M) {
  free(M->col);
  free(M->val);
  free(M);

  return(1);
}

/* convert from hash table to CSV matrix */
int csvPack(pCsv M) {
  csHash  *hm;
  double   val;
  int      i,j,idx;

  hm = M->hm;
  M->val = (double*)malloc((hm->nbe+1)*sizeof(double));
  M->col = (int*)malloc((hm->nbe+1)*sizeof(int));
  assert(M->val);
  assert(M->col);

  idx = M->nbe = 0;
  for (i=0; i<M->nr; i++) {
    for (j=0; j<M->nc; j++) {
      idx++;
      if ( cshGet(hm,i,j,&val) ) {
        M->val[idx] = val;
        M->col[idx] = idx;
        idx = 0;
        M->nbe++;
      }
    }
  }
  free(hm->item);
  free(M->hm);
  M->hm = 0;
  return(1);
}

/* add val to M(i,j) */
int csvPut(pCsv M,int i,int j,double val) {
  return(cshPut(M->hm,i,j,val));
}

/* set val to M(i,j) */
int csvSet(pCsv M,int i,int j,double val) {
  return(cshSet(M->hm,i,j,val));
}


/*======= Part 2: CSR matrix =======*/

/* alloc memory to store matrix */
pCsr csrNew(int nr,int nc,int nmax,char typ) {
  pCsr   M;

  M = (pCsr)calloc(1,sizeof(Csr));  
  assert(M);
  M->nr   = nr;
  M->nc   = nc;
  M->nbe  = 0;
  M->nmax = 0;
  M->hm   = cshNew(nr,nmax);
  M->typ  = typ;
  return(M);
}

/* alloc memory to store matrix */
void csrAlloc(pCsr M,int nr,int nc,int nmax,char typ) {
  M->nr   = nr;
  M->nc   = nc;
  M->nbe  = 0;
  M->nmax = 0;
  M->hm   = cshNew(nr,nmax);
  M->typ  = typ;
}

/* release memory */
int csrFree(pCsr M) {
  free(M->row);
  free(M->col);
  free(M->val);
  free(M);
  return(1);
}

/* optimize allocation for matrix structure */
static int csmPack(pCsr M) {
  if ( M->nbe < M->nmax ) {
    M->col = (int*)realloc(M->col,M->nbe*sizeof(int));
    M->val = (double*)realloc(M->val,M->nbe*sizeof(double));
    assert(M->col);
    assert(M->val);
    M->nmax = M->nbe;
  }
  return(1);
}

/* assign diagonal term M_ii to row[i] */
static void csmMii(pCsr M) {
  double   iv;
  int      i,j,ic;

  for (i=0; i<M->nr; i++) {
    if ( M->col[M->row[i]] == i )  continue;
    for (j=M->row[i]+1; j<M->row[i+1]; j++) {
      if ( M->col[j] == i ) {
        ic = M->col[M->row[i]];
        M->col[M->row[i]] = i;
        M->col[j] = ic;
        iv = M->val[M->row[i]];
        M->val[M->row[i]] = M->val[j];
        M->val[j] = iv;
        break;
      }
    }
  }
}

double *_val;
int    *_per,*_col;

static int compval(const void *a,const void *b) {
  int  diff = _col[*(int*)a] - _col[*(int*)b];
  return(0 < diff) - (diff < 0);
}

/* to improve access, sort M(i,j) ascending order */
static void csmSort(pCsr M){
  int    i,j,k,siz,nitm;

	/* memory alloc for sorting */
	siz = 64;
	_per = (int*)malloc(siz*sizeof(int));
	assert(_per);
	_col = (int*)malloc(siz*sizeof(int));
	assert(_col);
	_val = (double*)malloc(siz*sizeof(double));
	assert(_val);

  for (i=0; i<M->nr; i++) {
    nitm = M->row[i+1] - M->row[i];
		if ( nitm > siz ) {
			siz  = 1.5*nitm;
			_per = (int*)realloc(_per,siz*sizeof(int));
			assert(_per);
			_col = (int*)realloc(_col,siz*sizeof(int));
			assert(_col);
			_val = (double*)realloc(_val,siz*sizeof(double));
			assert(_val);
		}
		for (j=0; j<nitm; j++)  _per[j] = j;
		memcpy(_col,&M->col[M->row[i]],nitm*sizeof(int));
		memcpy(_val,&M->val[M->row[i]],nitm*sizeof(double));

    qsort(_per,nitm,sizeof(int),compval);

		for (j=M->row[i]; j<M->row[i+1]; j++) {
			k = j - M->row[i];
      M->col[j] = _col[_per[k]];
			M->val[j] = _val[_per[k]];
		}
  }
	free(_per);
	free(_col);
	free(_val);
}

/* convert from hash table to CSR matrix */
int csrPack(pCsr M) {
  csHash  *hm;
  hmat    *phm;
  int      i;

  if ( !M->hm )  return(csmPack(M));

  hm = M->hm;
  M->val = (double*)malloc(hm->nbe*sizeof(double));
  M->col = (int*)malloc(hm->nbe*sizeof(int));
  M->row = (int*)malloc((M->nr+1)*sizeof(int));
  assert(M->val);
  assert(M->col);
  assert(M->row);

  M->nbe = 0;
  assert(M->nr == hm->siz);
  for (i=0; i<M->nr; i++) {
    phm = &hm->item[i];
    if ( phm->j < 0 )  continue;
    M->row[i] = M->nbe;
    if ( fabs(phm->val) > CS_NUL ) {
      M->col[M->nbe] = phm->j;
      M->val[M->nbe] = phm->val;
      M->nbe++;
    }
    while ( phm->nxt ) {
      phm = &hm->item[phm->nxt];
      if ( fabs(phm->val) > CS_NUL ) {
        M->col[M->nbe] = phm->j;
        M->val[M->nbe] = phm->val;
        M->nbe++;
      }
    }
  }
  assert(M->nbe <= hm->nbe);
  M->row[M->nr] = M->nbe;
  free(hm->item);
  free(M->hm);
  M->hm = 0;
  M->nmax = M->nbe;

  /* optimization: sort ascending order of cols */
  csmSort(M);
  return(1);
}

/* add val to M(i,j) */
int csrPut(pCsr M,int i,int j,double val) {

  if ( fabs(val) < CS_NUL )  return(0);
  else if ( i >= M->nr || j >= M->nc )     return(0);
  else if ( (M->typ & CS_UT) && (j < i) )  return(0);
  else if ( (M->typ & CS_LT) && (j > i) )  return(0);

  if ( M->hm )
    return(cshPut(M->hm,i,j,val));
  else
    return(csmPut(M,i,j,val));
}

/* set val to M(i,j) */
int csrSet(pCsr M,int i,int j,double val) {
  if ( (M->typ & CS_UT) && (j < i) )  return(0);
  else if ( (M->typ & CS_LT) && (j > i) )  return(0);

  if ( M->hm )
    return(cshSet(M->hm,i,j,val));
  else
    return(csmSet(M,i,j,val));
}


/*======= Part 3: multi-threaded matrix operators =======*/

/* initialize parallel lib with ncpu procs */
void csrInit(int ncpu) {

  CSR_libCpu = ncpu;
  if ( ncpu > CS_MAXCPU ) {
    fprintf(stdout,"  ## Number of Cpu's limited to %d\n",CS_MAXCPU);
    CSR_libCpu = CS_MAXCPU;
  }
  CSR_libId = InitParallel(ncpu);
  assert(CSR_libId);
}

void csrStop() {
  StopParallel(CSR_libId);
  CSR_libId  = 0;
  CSR_libCpu = 1;
}


/* matrix transpose Mt = A^t */
pCsr csrTr(pCsr M) {
  pCsr    Mt;
  int    *wt,i,j,ic,nbe;

  /* memory allocation */
  Mt = csmNew(M->nc,M->nr,M->nbe,0);
  if ( M->typ & CS_UT )  Mt->typ |= CS_LT;
  else if ( M->typ & CS_LT )  Mt->typ |= CS_UT;
  else if ( M->typ & CS_SYM ) Mt->typ |= CS_SYM;

  /* count number of entries / columns in M */
  wt = (int*)calloc(M->nc,sizeof(int));
  assert(wt);
  for (j=0; j<M->nbe; j++) 
   if ( fabs(M->val[j]) >= CS_NUL )  wt[M->col[j]]++;

  /* set rows pointers in Mt */
  nbe = 0;
  for (i=0; i<Mt->nr; i++) {
    Mt->row[i] = nbe;
    nbe  += wt[i];
    wt[i] = Mt->row[i];
  }
  Mt->row[Mt->nr] = nbe;

  /* store Mt(i,j) = M(j,i) */
  for (j=0; j<M->nr; j++) {
    for (i=M->row[j]; i<M->row[j+1]; i++) {
      if ( fabs(M->val[i]) < CS_NUL )  continue;
      ic = wt[M->col[i]]++;   /* next entry */
      Mt->col[ic] = j;        /* Mt(j,i) = M(i,j) */
      Mt->val[ic] = M->val[i];
    }
  }
  free(wt);
  return(Mt);
}

/* C = l*A + m*B */
pCsr csrAdd(pCsr A,pCsr B,double l,double m) {
  pCsr     C;
  double  *x;
  int     *wt,i,j,nr,nc,ic,nbe;

  if ( A->hm || B->hm )  return(0);

  /* resulting matrix */
  nr = CS_MIN(A->nr,B->nr);
  nc = CS_MIN(A->nc,B->nc);
  C  = csmNew(nr,nc,A->nbe + B->nbe,0);
  if ( (A->typ & CS_SYM) && (B->typ & CS_SYM) )  C->typ |= CS_SYM;
  if ( (A->typ & CS_UT)  && (B->typ & CS_UT) )   C->typ |= CS_UT;
  if ( (A->typ & CS_LT)  && (B->typ & CS_LT) )   C->typ |= CS_LT;

  /* memory alloc */
  wt = (int*)calloc(nc,sizeof(int));
  assert(wt);
  x  = (double*)calloc(nc,sizeof(double));
  assert(x);

  nbe = 0;
  for (i=0; i<nr; i++) {
    C->row[i] = nbe;
    /* process A(i,:) */
    for (j=A->row[i]; j<A->row[i+1]; j++) {
      ic = A->col[j];
      if ( fabs(A->val[j]) < CS_NUL )  continue;
      else if ( wt[ic] == 0 ) {
        C->col[nbe++] = ic;
        wt[ic]++;
        x[ic] = l * A->val[j];
      }
      else {
        x[ic] += l * A->val[j];
      }
    }
    /* process B(i,:) */
    for (j=B->row[i]; j<B->row[i+1]; j++) {
      ic = B->col[j];
      if ( fabs(B->val[j]) < CS_NUL )  continue;
      else if ( wt[ic] == 0 ) {
        C->col[nbe++] = ic;
        wt[ic]++;
        x[ic] = m * B->val[j];
      }
      else {
        x[ic] += m * B->val[j];
      }
    }
    /* store in C(i,:) */
    for (j=C->row[i]; j<nbe; j++) {
      C->val[j] = x[C->col[j]];
    }
    memset(wt,0,nc*sizeof(int));
    memset(x,0,nc*sizeof(double));
  }
  free(x);
  free(wt);
  C->row[C->nr] = nbe;
  C->nbe = nbe;
  C->col = realloc(C->col,nbe*sizeof(int));
  C->val = realloc(C->val,nbe*sizeof(double));

  if ( C->typ & CS_SYM )  csmMii(C);
  return(C);
}

/* C = A x At, store result in symmetric matrix */
pCsr csrMulAAt(pCsr A) {
  pCsr     C;
  double  *x,cij;
  int      i,j,k,nbe,nbemax;

  if ( A->hm )  return(0);

  /* resulting matrix is symmetric */
  nbemax = 20*A->nbe;
  C = csmNew(A->nr,A->nr,nbemax,CS_UT+CS_SYM);

  /* memory alloc */
  x = (double*)calloc(A->nc,sizeof(double));
  assert(x);

  nbe = 0;
  for (i=0; i<A->nr; i++) {
    /* store A(i,:) sorted column values for row i */
    for (j=A->row[i]; j<A->row[i+1]; j++)
      x[A->col[j]] = A->val[j];

    C->row[i] = nbe;
    /* process At(:,j) = A(j,:) */
    for (j=i; j<A->nr; j++) {
      cij = 0.0;
      for (k=A->row[j]; k<A->row[j+1]; k++) {
        cij += x[A->col[k]] * A->val[k];
      }
      if ( fabs(cij) < CS_NUL )  continue;
      C->col[nbe] = j;
      C->val[nbe] = cij;
      nbe++;
      /* reallocation */
      if ( nbe >= nbemax ) {
        nbemax = (int)(1.5*nbemax);
        C->col = (int*)realloc(C->col,nbemax*sizeof(int));
        assert(C->col);
        C->val = (double*)realloc(C->val,nbemax*sizeof(double));
        assert(C->val);
      }
    }
    /* reset vector columns */
    memset(x,0,A->nc*sizeof(double));
  }
  free(x);
  C->row[C->nr] = nbe;
  C->nbe = nbe;
  C->col = (int*)realloc(C->col,nbe*sizeof(int));
  C->val = (double*)realloc(C->val,nbe*sizeof(double));

  csmMii(C);
  return(C);
}

/* C(m,p) = A(m,n) x B(n,p) (non symmetric matrices) */
pCsr csrMul(pCsr A,pCsr B) {
  pCsr     C;
  double  *x,v;
  int     *wt,i,j,k,l,ic,nbe;

  if ( A->hm || B->hm )  return(0);
  assert((A->nc==B->nr));

  /* resulting matrix */
  C = csmNew(A->nr,B->nc,A->nbe + B->nbe,0);
  wt = (int*)malloc(A->nr*sizeof(int));
  assert(wt);
  x  = (double*)malloc(A->nr*sizeof(double));
  assert(x);

  nbe = 0;
  for (j=0; j<B->nc; j++) {
    C->row[j] = nbe;
    for (i=B->row[j]; i<B->row[j+1]; i++) {
      ic = B->col[i];
      v  = B->val[i];
      for (k=A->row[ic]; k<A->row[ic+1]; k++) {
        l = A->col[k];
        if ( wt[l] < j+1 ) {
          wt[l] = j+1;
          nbe++;
          C->col[nbe] = l;
          x[l] = v * A->val[k];
        }
        else
          x[l] += v * A->val[k];
      }
    }
    for (i=C->row[j]; i<nbe; i++)
      C->val[i] = x[C->col[i]];
  }
  C->row[B->nc] = nbe;

  free(wt);
  free(x);
  C->nbe = nbe;
  C->col = realloc(C->col,nbe*sizeof(int));
  C->val = realloc(C->val,nbe*sizeof(double));

  if ( C->typ & CS_SYM )  csmMii(C);
  return(C);
}

/*--------------- multithreaded part: local ---------------*/
/* compute z = l.Ax + m.y */
static void csr_axpy(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  pCsr    A;
  double *x,*y,*z,li;
  int     i,j,ic;

  A = arg->A;
  x = arg->x;
  y = arg->y;
  z = arg->z;
  for (i=startAdr-1; i<stopAdr; i++) {
    li = 0.0;
    for (j=A->row[i]; j<A->row[i+1]; j++)
      li += A->val[j] * x[A->col[j]];
    z[i] = arg->l*li + arg->m*y[i];
  }
  if ( A->typ & CS_SYM ) {
    /* use optim: M_ii in M->row[i] */
    for (i=startAdr-1; i<stopAdr; i++) {
      for (j=A->row[i]+1; j<A->row[i+1]; j++) {
        ic = A->col[j];
				li = arg->l*(A->val[j] * x[i]);
        z[ic] += li;
      }
    }
  }
}

/* compute y = A.x */
static void csr_ax(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  pCsr    A;
  double *x,*y,dd;
  int     i,j,ic;

  A = arg->A;
  x = arg->x;
  y = arg->y;
  for (i=startAdr-1; i<stopAdr; i++) {
    y[i] = 0.0;
    for (j=A->row[i]; j<A->row[i+1]; j++) {
			dd    = A->val[j] * x[A->col[j]];
      y[i] += dd;
	  }
  }
  if ( A->typ & CS_SYM ) {
    /* use optim: M_ii in M->row[i] */
    for (i=startAdr-1; i<stopAdr; i++) {
      for (j=A->row[i]+1; j<A->row[i+1]; j++) {
        ic = A->col[j];
				dd     = A->val[j] * x[i];
        y[ic] += dd;
      }
    }
  }
}

/* compute y = A^t.x */
static void csr_atx(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  pCsr    A;
  double *x,*y,dd;
  int     i,j,ic;

  A = arg->A;
  x = arg->x;
  y = arg->y;
  for (j=startAdr-1; j<stopAdr; j++) {
    for (i=A->row[j]; i<A->row[j+1]; i++) {
      ic = A->col[i];
      y[ic] = y[ic] + A->val[i] * x[j];
    }
  }
  if ( A->typ & CS_SYM ) {
    /* use optim: M_ii in M->row[i] */
    for (j=startAdr-1; j<stopAdr; j++) {
      dd = 0.0;
      for (i=A->row[j]+1; i<A->row[j+1]; i++)
        //y[j] = y[j] + A->val[i] * x[A->col[i]];
        dd += A->val[i] * x[A->col[i]];
      y[j] += dd;
    }
  }
}

static void csr_lxmy(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  double  *x,*y,*z,l,m,dd;
  int     i;

  x = arg->x;
  y = arg->y;
	z = arg->z;
  l = arg->l;
	m = arg->m;
  for (i=startAdr-1; i<stopAdr; i++) {
    dd   = l*x[i] + m*y[i];
		z[i] = dd;
	}
}

static void csr_lxy(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  double  *x,*y,l;
  int     i;

  x = arg->x;
  y = arg->y;
  l = arg->l;
  for (i=startAdr-1; i<stopAdr; i++)
    y[i] = l * x[i];
}

static void csr_xy(int startAdr,int stopAdr,int PthIdx,CsrArg *arg) {
  double  r,*x,*y;
  int     i;

  x = arg->x;
  y = arg->y;
  r = 0.0;
  for (i=startAdr-1; i<stopAdr; i++)
    r += x[i] * y[i];
  arg->r[PthIdx] = r;
}

/*--------------- multithreaded part: global ---------------*/
/* multithreaded: y = A*x */
int csrAx(pCsr A,double *x,double *y) {
  CsrArg   arg;
  float    acc;
  int      typid;

  if ( !x || !y )  return (0);
  arg.A = A;
  arg.x = x;
  arg.y = y;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,A->nr);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_ax,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_ax(1,A->nr,0,&arg);
  return(1);
}

/* multithreaded: y = A^t*x (A transposed directly) */
int csrAtx(pCsr A,double *x,double *y) {
  CsrArg   arg;
  float    acc;
  int      typid;

  if ( !x || !y )  return (0);
  memset(y,0,A->nc*sizeof(double));
  arg.A = A;
  arg.x = x;
  arg.y = y;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,A->nr);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_atx,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_atx(1,A->nr,0,&arg);
  return(1);
}

/* multithreaded: z = l.A*x + m.y */
int csrAxpy(pCsr A,double *x,double *y,double *z,double l,double m) {
  CsrArg   arg;
  float    acc;
  int      typid;

  if ( !x || !y || !z )  return (0);
  arg.A = A;
  arg.x = x;  arg.y = y;  arg.z = z;
  arg.l = l;  arg.m = m;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,A->nr);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_axpy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_axpy(1,A->nr,0,&arg);
  return(1);
}

/* multithreaded: z = l.A^t*x + m.y */
int csrAtxpy(pCsr A,double *x,double *y,double *z,double l,double m) {
  CsrArg   arg;
  float    acc;
  int      typid;

  if ( !x || !y || !z )  return (0);
  arg.A = A;
  arg.x = x;
  arg.y = z;
  if ( CSR_libId ) {
		/* z = arg.y = A^t*x */
    typid = NewType(CSR_libId,A->nr);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_atx,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);

    /* z = l.z + m.y = l.A^t*x + m.y */
    arg.x = z;  arg.y = y;  arg.z = z;
		arg.l = l;	arg.m = m;
    typid = NewType(CSR_libId,A->nc);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_lxmy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else {
		/* z = arg.y = A^t*x */
    csr_atx(1,A->nr,0,&arg);
    arg.x = z;  arg.y = y;  arg.z = z;
		arg.l = l;  arg.m = m;
    /* z = l.z + m.y */
    csr_lxmy(1,A->nc,0,&arg);
  }
  
  //for (i=0; i<A->nc; i++)  z[i] = l*z[i] + m*y[i];
  return(1);
}

/* res = <y=A*x,x> */
double csrAxdotx(pCsr A,double *x,double *y) {
  CsrArg   arg;
  double   axx;
  float    acc;
  int      i,typid;

  if ( !x || !y )  return (0);
  memset(y,0,A->nr*sizeof(double));
  arg.A = A;
  arg.x = x;
  arg.y = y;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,A->nr);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_ax,(void *)&arg);
    assert(acc);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_xy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else {
    csr_ax(1,A->nr,0,&arg);
    csr_xy(1,A->nr,0,&arg);
  }
  axx = 0.0;
  for (i=0; i<CSR_libCpu; i++)  axx += arg.r[i];
  return(axx);
}

/* z[i] = l*x[i] + m*y[i] */
void csrlXmY(double *x,double *y,double *z,double l,double m,int n) {
  CsrArg   arg;
	float    acc;
  int      typid;

  arg.x = x;
  arg.y = y;
	arg.z = z;
	arg.l = l;
	arg.m = m;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,n);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_lxmy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_lxmy(1,n,0,&arg);
}

/* y[i] = l*x[i] */
void csrlX(double *x,double *y,double l,int n) {
  CsrArg   arg;
	float    acc;
  int      typid;

  arg.x = x;
  arg.y = y;
	arg.l = l;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,n);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_lxy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_lxy(1,n,0,&arg);
}

/* r = <x,y> */
double csrXY(double *x,double *y,int n) {
  CsrArg    arg;
  double    xy;
  float     acc;
  int       i,typid;

  arg.x = x;
  arg.y = y;
  if ( CSR_libId ) {
    typid = NewType(CSR_libId,n);
    assert(typid);
    acc = LaunchParallel(CSR_libId,typid,0,(void *)csr_xy,(void *)&arg);
    assert(acc);
    FreeType(CSR_libId,typid);
  }
  else
    csr_xy(1,n,0,&arg);

  xy = 0.0;
  for (i=0; i<CSR_libCpu; i++)  xy += arg.r[i];
  return(xy);
}


/* 1-norm of CSR matrix = max (sum(abs(A))) */
double csrNorm(pCsr M) {
  int      i,j;
  double   nn,s;
  
  if ( !CS_CSR(M) || !M->val )  return(-1);

  nn = 0.0;
  for (i=0; i<M->nr; i++) {
    s = 0.;
    for (j=M->row[i]; j<M->row[i+1]; j++) 
      s += fabs(M->val[j]);
    nn = CS_MAX(nn,s);
  }

  return (nn);
}

int csrPrint(pCsr A,int brief) {
  int    i,j;

  if ( !A ) { 
    fprintf(stdout,"(null)\n");
    return(0);
  }
  if ( brief < 0 ) {
    fprintf(stdout,"  Sparse Matrix Library, Version %s (%s)\n",CS_VER,CS_REL);
    fprintf(stdout,"  %s\n",CS_CPY);
  }
  if ( A->nbe < 0 ) {
    fprintf(stdout,"  %d-by-%d, nbe: %d, 1-norm: %g\n",A->nr,A->nc,A->nbe,csrNorm(A));
    for (i=0; i<A->nr; i++) {
      fprintf(stdout,"    row %d: loc %d to %d\n",i,A->row[i],A->row[i+1]-1);
      for (j=A->row[i]; j<A->row[i+1]; j++) {
        fprintf(stdout,"     %4d : %E\n", A->col[j],A->val ? A->val[j] : 1);
        if ( brief && i>20 ) {
          fprintf(stdout,"  ...\n"); 
          return (1);
        }
      }
    }
  }
  else {
    fprintf(stdout,"\n  raw: %d-by-%d, nbe: %d\n",A->nr,A->nc,A->nbe);
    for (i=0; i<A->nr; i++) {
      fprintf(stdout,"    row [%d]\n      ",i);
      for (j=A->row[i]; j<A->row[i+1]; j++) {
        fprintf(stdout," %d: %g  ",A->col[j],A->val[j]);
        fflush(stdout);
      }
      fprintf(stdout,"\n");
      if ( brief && i > 9 )  break;
    }
  }

  return (1);
}

void csrPrLine(pCsr A,int i) {
  int   k,l;

  fprintf(stdout,"row[%d]:\n",i);
  l = 0;
  for (k=A->row[i]; k<A->row[i+1]; k++) {
    fprintf(stdout,"  %6d: %+e",A->col[k],A->val[k]);
    fflush(stdout);
    ++l;
    if ( l % 5 == 0 ) {
      fprintf(stdout,"\n");
      l = 0;
    }
  }
  fprintf(stdout,"\n");
}

void csrPrVal(pCsr A,int i,int j) {
  int   k;

  for (k=A->row[i]; k<A->row[i+1]; k++) {
    if ( A->col[k] == j ) {
      fprintf(stdout,"A[%d][%d] = %f\n",i,j,A->val[k]);
      return;    
    }
  }
}

#define NBL   10

pCsr csrLoad(char *name) {
  FILE  *in;
	pCsr   A;
  int    i;

  in = fopen(name,"r");
  assert(in);

	A = (pCsr)calloc(1,sizeof(Csr));
	fscanf(in,"%d %d %d",&A->nr,&A->nc,&A->nbe);
	fscanf(in,"%c",&A->typ);

	A->val = (double*)malloc(A->nbe*sizeof(double));
	A->col = (int*)malloc(A->nbe*sizeof(int));
	A->row = (int*)malloc((A->nr+1)*sizeof(int));
	assert(A->val);
	assert(A->col);
	assert(A->row);

  for (i=0; i<=A->nr; i++) {
    fscanf(in,"%d",&A->row[i]);
  }

  for (i=0; i<A->nbe; i++) {
    fscanf(in,"%d",&A->col[i]);
	}

  for (i=0; i<A->nbe; i++) {
    fscanf(in,"%lg",&A->val[i]);
  }

	fclose(in);
	return(A);
}


void csrSave(pCsr A,char *name) {
  FILE  *out;
  int    i,l;

  out = fopen(name,"w");
  assert(out);

  fprintf(out,"%d %d %d\n",A->nr,A->nc,A->nbe);
	fprintf(out,"%c\n",A->typ);
  l = 0;
	for (i=0; i<=A->nr; i++) {
		fprintf(out,"%d ",A->row[i]);
	  if ( ++l == NBL ) {
	  	fprintf(out,"\n");
			l = 0;
	  }
  }
  if ( l )  fprintf(out,"\n");

	l = 0;
  for (i=0; i<A->nbe; i++) {
    fprintf(out,"%d ",A->col[i]);
	  if ( ++l == NBL ) {
	    fprintf(out,"\n");
		  l = 0;
	  }
  }
  if ( l )  fprintf(out,"\n");

  l = 0;
	for (i=0; i<A->nbe; i++) {
    fprintf(out,"%g ",A->val[i]);
	  if ( ++l == NBL ) {
		  fprintf(out,"\n");
		  l = 0;
	  }
	}
  if ( l )  fprintf(out,"\n");

  fclose(out);
}

