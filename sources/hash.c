#include "elastic.h"


/* hashing structure */
typedef struct {
  int   min,max,nod,nxt;
} Cell;

typedef struct {
  Cell   *cell;
  int     hsiz,nmax,hnxt;
} Htab;


/* return P2 node along edge (a,b) if created */
static int nodeP2(Htab *ht,int a,int b) {
  Cell   *pc;
  int     sum,min,max;

  sum = (a+b) % ht->hsiz;
  pc  = &ht->cell[sum];

  if ( !pc->min )  return(0);

  /* check linked list */
  min = LS_MIN(a,b);
  max = LS_MAX(a,b);
  do {
    if ( pc->min == min && pc->max == max )  return(pc->nod);
    else if ( !pc->nxt )  return(0);
    pc = &ht->cell[pc->nxt];
  }
  while ( pc );

  return(0);
}


/* create new P2 node on (a,b) or return index */
static int hedge(Htab *ht,int a,int b,int *na) {
  Cell   *pc;
	int     j,min,max,sum;

  sum = (a+b) % ht->hsiz;
  pc  = &ht->cell[sum];
  min = LS_MIN(a,b);
  max = LS_MAX(a,b);

  /* insert new midpoint */
	if ( !pc->min ) {
    *na = *na + 1;
    pc->min = min;
    pc->max = max;
    pc->nod = *na;
    pc->nxt= 0;
    return(*na);
  }
  /* check if existing */
  else if ( pc->min == min && pc->max == max )
    return(pc->nod);
  else {
    while ( (pc->nxt > 0) && (pc->nxt < ht->nmax) ) {
      pc = &ht->cell[pc->nxt];
			if ( (pc->min == min) && (pc->max == max) )  return(pc->nod);
    }
    /* insert new midpoint */
    pc->nxt = ht->hnxt;
    pc      = &ht->cell[pc->nxt];
    ++ht->hnxt;

    *na = *na + 1;
    pc->min = min;
    pc->max = max;
    pc->nod = *na;
    pc->nxt = 0;

    /* check for overflow */
    if ( ht->hnxt == ht->nmax ) {
	    ht->nmax = 1.2 * ht->nmax;
	    ht->cell = (Cell*)realloc(ht->cell,ht->nmax*sizeof(Cell));
      assert(ht->cell);
      for (j=ht->hnxt; j<ht->nmax; j++)  ht->cell[j].nxt = j+1;
    }
		return(pc->nod);
  }
  return(0);
}


/* hash mesh edges for creating P2 nodes */
int hashar_3d(LSst *lsst) {
	Htab    ht;
  pTetra  pt;
  int     k,na,ip;
	char    i,i1,i2,off;
  static int edg[6][2] = {0,1, 0,2, 0,3, 1,2, 1,3, 2,3};

  /* alloc hash table */
  ht.nmax = (int)(8.2*lsst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = lsst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;

  /* loop over tetrahedra */
	na  = 0;
  off = 4;
	for (k=1; k<=lsst->info.ne; k++) {
		pt = &lsst->mesh.tetra[k];
    /* loop over edges */
		for (i=0; i<6; i++) {
			i1 = edg[i][0];
			i2 = edg[i][1];
			ip = hedge(&ht,pt->v[i1],pt->v[i2],&na);
      pt->v[i+off] = ip;
		}
	}
  lsst->info.na = na;
  free(ht.cell);

  return(na);
}


/* hash mesh edges for creating P2 nodes */
int hashar_2d(LSst *lsst) {
  Htab    ht;
  pTria   pt;
  int     k,na,ip;
	char    i,i1,i2,off;

  /* alloc hash table */
  ht.nmax = (int)(3.2*lsst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = lsst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;

  /* loop over triangles */
  na  = 0;
  off = 3;
  for (k=1; k<=lsst->info.nt; k++) {
    pt = &lsst->mesh.tria[k];
    for (i=0; i<3; i++) {
      i1 = (i+1) % 3;
      i2 = (i+2) % 3;
      ip = hedge(&ht,pt->v[i1],pt->v[i2],&na);
      pt->v[i+off] = ip;
    }
  }
  lsst->info.na = na;
  free(ht.cell);

	return(na);
}