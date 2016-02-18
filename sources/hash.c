#include "elastic.h"

#define KTA     7
#define KTB    11


/* get P2 node */
int hashP2(Hash *hash,int a,int b) {
  hedge  *ph;
  int     key,ia,ib;

  ia  = LS_MIN(a,b);
  ib  = LS_MAX(a,b);
  key = (KTA*ia + KTB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->ia )  return(0);
  if ( ph->ia == ia && ph->ib == ib )  return(ph->k);
  while ( ph->nxt ) {
    ph = &hash->item[ph->nxt];
    if ( ph->ia == ia && ph->ib == ib )  return(ph->k);
  }

  return(0);
}


static int hashPut(Hash *hash,int a,int b,int *na) {
  hedge  *ph;
	int     j,ia,ib,key;

	ia = LS_MIN(a,b);
	ib = LS_MAX(a,b);
  key = (KTA*ia + KTB*ib) % hash->siz;
  ph  = &hash->item[key];

	if ( ph->ia == ia && ph->ib == ib )
		return(0);
	else if ( ph->ia ) {
    while ( ph->nxt && ph->nxt < hash->max ) {
      ph = &hash->item[ph->nxt];
			if ( ph->ia == ia && ph->ib == ib )  return(0);
    }
    /* insert new midpoint */
		*na = *na + 1;
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
		ph->ia  = ia;
		ph->ib  = ib;
		ph->k   = *na;
		ph->nxt = 0;
    ++hash->nxt;
    /* check for overflow */
    if ( hash->nxt >= hash->max ) {
	    hash->max *= 1.2;
	    hash->item  = (hedge*)realloc(hash->item,hash->max*sizeof(hedge));
      assert(hash->item);
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
		return(1);
  }
  /* insert new midpoint */
	*na = *na + 1;
  ph->ia = ia;
  ph->ib = ib;
  ph->k  = *na;
  ph->nxt= 0;

	return(1);
}


/* hash mesh edges for creating P2 nodes */
int hashar_3d(LSst *lsst) {
	pTetra  pt;
  int     k;
	char    i,i1,i2;
  static int edg[6][2] = {0,1, 0,2, 0,3, 1,2, 1,3, 2,3};

  /* adjust hash table params */
  lsst->hash.siz  = lsst->info.np;
  lsst->hash.nxt  = lsst->info.np;
  lsst->hash.max  = 3.2*lsst->info.np + 1;
  lsst->hash.item = (hedge*)calloc(lsst->hash.max,sizeof(hedge));
  assert(lsst->hash.item);
  for (k=lsst->info.np; k<lsst->hash.max; k++)  lsst->hash.item[k].nxt = k+1;

  /* loop over tetrahedra */
	lsst->info.na = 0;
	for (k=1; k<=lsst->info.ne; k++) {
		pt = &lsst->mesh.tetra[k];
    /* loop over edges */
		for (i=0; i<6; i++) {
			i1 = edg[i][0];
			i2 = edg[i][1];
			hashPut(&lsst->hash,pt->v[i1],pt->v[i2],&lsst->info.na);
		}
	}

	return(lsst->info.na);
}


/* hash mesh edges for creating P2 nodes */
int hashar_2d(LSst *lsst) {
	pTria   pt;
  int     k;
	char    i,i1,i2;

  /* adjust hash table params */
  lsst->hash.siz  = lsst->info.np;
  lsst->hash.nxt  = lsst->info.np;
  lsst->hash.max  = 3.2*lsst->info.np + 1;
  lsst->hash.item = (hedge*)calloc(lsst->hash.max,sizeof(hedge));
  assert(lsst->hash.item);
  for (k=lsst->info.np; k<lsst->hash.max; k++)  lsst->hash.item[k].nxt = k+1;

  /* loop over triangles */
	lsst->info.na = 0;
	for (k=1; k<=lsst->info.nt; k++) {
		pt = &lsst->mesh.tria[k];
		for (i=0; i<3; i++) {
			i1 = (i+1) % 3;
			i2 = (i+2) % 3;
			hashPut(&lsst->hash,pt->v[i1],pt->v[i2],&lsst->info.na);
		}
	}

	return(lsst->info.na);
}