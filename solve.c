#include "sparse.h"


#define CS_EPSD    1.0e-30
#define CS_EPSD2   1.0e-200
#define CS_KRYL    100
#define CS_KRY1    CS_KRYL*(CS_KRYL+3) / 2
#define CS_OMEGA   1.2   /* >1 to speed up convergence of slow-converging process */


/* SSOR preconditioner to solve (L*d*L^t)x = b
 * input:  L=A^t lower triangular, x,b  vectors
 * output: x modified  */
int csrSSOR(pCsr A,pCsr L,double *x,double *b) {
  double  *dia,sum,dd,omega;
  int      i,j;

  omega = CS_OMEGA;
  dia = (double*)calloc(L->nr,sizeof(double));
  assert(dia);

  for (i=0; i<L->nr; i++) {
    sum = 0.0;
    for (j=L->row[i]; j<L->row[i+1]-1; j++)
      sum += L->val[j] * x[L->col[j]];
    dia[i] = L->val[j];     /* optim: diag term is in last position */
		if ( fabs(dia[i]) > CS_EPSD2 )
      x[i] = (b[i] - sum) * omega / dia[i];
  }

  /* diagonal terms */
  dd  = (2.0 - omega) / omega;
  for (i=0; i<L->nr; i++)
    x[i] = x[i] * dd * dia[i];

  for (i=A->nr-1; i>=0; i--) {
    sum = 0.0;
    for (j=A->row[i]+1; j<A->row[i+1]; j++)
       sum += A->val[j] * x[A->col[j]];
    if ( fabs(dia[i]) > CS_EPSD2 )
		  x[i] = (x[i] - sum) * omega / dia[i];
  }
  free(dia);
  return(1);
}

/* solve Ax=b, using gradient. Return code:
   >0 upon completion, <0 no convergence
   0 memory problem 
  -1 ill conditioned matrix
  -2 max it reached  */
int csrGradient(pCsr A,double *x,double *b,double *er,int *ni) {
  double   *y,*ay,dp,rm,rmp,tol,alpha,err;
  int       i,n,it,ier,nit;

  if ( !x || !b )  return(0);
  n = A->nr;
  y = (double*)calloc(n,sizeof(double));
  if ( !y )  return(0);

  /* compute r0 = b - A x0 */
  memcpy(y,b,n*sizeof(double));
  if ( !csrAxpy(A,x,y,NULL,-1.,1.) ) {
    free(y);
    return(0);
  }

  err = *er;
  nit = *ni;
  it  = 0;
  ier = 1;
  rm  = 0.0;
  for (i=0; i<n; i++)  rm += y[i] * y[i];
  rmp = rm;
  if ( fabs(rmp) < CS_EPSD2 ) {
    free(y);
    return(1);
  }
  else if ( rmp > CS_TGV2 )
    rmp /= CS_TGV2;
  tol = err*err*rm;

  ay = (double*)malloc(n*sizeof(double));
  if ( !ay ) {
    free(y);
    return(0);
  }

  if ( nit < 0 )  nit = CS_MAXIT;
  while ( (tol < rm) && (it < nit) ) {
    dp = csrAxdotx(A,y,ay);
    if ( fabs(dp) <= CS_EPSD2 )  break;

    /* x = x-alpha.y; y = y-alpha.<A,y> */
    alpha = rm / dp;
    rm = 0.0;
    for (i=0; i<n; i++) {
      x[i] -= alpha * y[i];
    }
    for (i=0; i<n; i++) {
      y[i] -= alpha * ay[i];
      rm   += y[i]*y[i];
    }

    if ( rm > 2.0*rmp ) {
      ier = -1;
      break;
    }
    rmp = rm;
    it++;
  }
  if ( it > nit )  ier = -2;
  *er = sqrt(rm);
  *ni = it;

  free(y);
  free(ay);
  return(ier);
}


/* solve Ax=b, using conjugate gradient. Return code:
   >0 upon completion, <0 no convergence
   0 memory problem 
  -1 ill conditioned matrix
  -2 max it reached  */
int csrConjGrad(pCsr A,double *x,double *b,double *er,int *ni) {
  double   *y,*ap,*p,dp,nn,rm,rmp,rmn,alpha,beta,err;
  int       i,n,it,ier,nit;
  
  if ( !x || !b )  return(0);
  n = A->nr;
  y = (double*)malloc(n*sizeof(double));
  assert(y);
  nn = csrXY(x,x,A->nr);

  /* compute R0 = b - A.x0 */
  if ( nn < CS_EPSD2 ) {
    memcpy(y,b,A->nr*sizeof(double));
  }
  else {
    csrAxpy(A,x,b,y,-1.,1.);
  }
  rm  = csrXY(y,y,n);
  rmp = rm;
  if ( fabs(rmp) < CS_EPSD2 ) {
    free(y);
    return(1);
  }
  else if ( rmp > CS_TGV2 )
    rmp /= CS_TGV2;

  p  = (double*)malloc(n*sizeof(double));
  ap = (double*)malloc(n*sizeof(double));
  assert(p);
  assert(ap);

  /* p1 = r0 */
  memcpy(p,y,n*sizeof(double));
  err = *er;
  err = err*err*rmp;
  nit = *ni;
  it  = 0;
  ier = 1;
  while ( (err < rm) && (++it <= nit) ) {
    /* alpha_m = <R_m-1,R_m-1> / <AP_m,P_m> */
    dp = csrAxdotx(A,p,ap);
    if ( fabs(dp) <= CS_EPSD2 )  break;

    /* X_m = X_m-1 + alpha_m.P_m; R_m = R_m-1 - alpha_m AP_m */ 
    alpha = rm / dp;
    rmn   = 0.0;
		csrlXmY(x,p,x,1.0,alpha,n);
		csrlXmY(y,ap,y,1.0,-alpha,n);
		rmn = csrXY(y,y,n);
		/*
    for (i=0; i<n; i++) {
      x[i] = x[i] + alpha * p[i];
      y[i] = y[i] - alpha * ap[i];
      rmn  += y[i] * y[i];
    }
    */
    /* P_m+1 = R_m + beta_m P_m */
    beta = rmn / rm;
    for (i=0; i<n; i++) {
      p[i] = y[i] + beta * p[i];
    }
    rm = rmn;
    // printf("  GC1: it  %d err %E   rm %E\n",it,err,rmn);
  }
  if ( it > nit )  ier = -2;
  *er = sqrt(rm / rmp);
  *ni = it;

  free(y);
  free(p);
  free(ap);
  return(ier);
}


/* GC modified */
int csrConjGradGen(pCsr A,double *x,double *b,double *ud,char *udt,char nohom,double *er,int *ni) {
  double   *y,*ap,*p,dp,rm,rmp,rmn,alpha,beta,err;
  int       i,n,it,ier,nit;

  if ( !x || !b )  return(0);
  if ( nohom && (!ud || !udt ) )   return(0);

  n = A->nr;
  y = (double*)malloc(n*sizeof(double));
  if ( !y )  return(0);

  memcpy(y,b,n*sizeof(double));
  if ( nohom )  csrAxpy(A,ud,b,y,-1.0,1.0);
  
  /* modify x0 */
  for (i=0; i<n; i++) {
    if ( udt[i] )  x[i] = 0.0; 
  }

  /* compute R0 = y - A.x0 */
  csrAxpy(A,x,y,y,-1.,1.);
  /* modify y, compute R0*R0 */  
  rmp = 0.0;
  for (i=0; i<n; i++) {
    if ( udt[i] )  
      y[i] = 0.0;
    else 
      rmp += y[i] * y[i];
  }

  if ( rmp < CS_EPSD ) {
    if ( nohom ) 
      for (i=0; i<n; i++)  if ( udt[i] ) x[i] = ud[i];
    *er = rmp;
    *ni = 0;
    free(y);
    return(1);
  }
  else if ( rmp > CS_TGV2 )
    rmp /= CS_TGV2;

  p =(double*)malloc(n*sizeof(double));
  assert(p);
  ap = (double*)malloc(n*sizeof(double));
  assert(ap);

  /* p1 = r0 */
  memcpy(p,y,n*sizeof(double));
  it  = 0;
  ier = 1;
  err = *er;
  err = err*err*rmp;
  nit = *ni;
  rm  = rmp;
  if ( nit < 0 )  nit = CS_MAXIT;
  while ( (err < rm) && (++it <= nit) ) {
    /* alpha_m = <R_m-1,R_m-1> / <AP_m,P_m> */
    dp = csrAxdotx(A,p,ap);
    if ( fabs(dp) <= CS_EPSD2 )  break;

    /* X_m = X_m-1 + alpha_m.P_m; R_m = R_m-1 - alpha_m AP_m */ 
    alpha = rm / dp;
    rmn   = 0.0;
    for (i=0; i<n; i++) {
      if ( udt[i] )  ap[i] = 0.0;
      x[i] += alpha * p[i];
      y[i] -= alpha * ap[i];
      rmn  += y[i] * y[i];
    }

    /* P_m+1 = R_m + beta_m P_m */
    if ( fabs(rm) <= CS_EPSD )  break;
    beta = rmn / rm;
    for (i=0; i<n; i++)
      p[i] = y[i] + beta * p[i];
    rm = rmn;
  }
  if ( nohom ) 
    for (i=0; i<n; i++)  if ( udt[i] )  x[i] = ud[i];

  if ( it > nit )  ier = -2;
  *er = sqrt(rm / rmp);
  *ni = it;

  free(y);
  free(p);
  free(ap);
  return(ier);
}


/* solve Ax = b using preconditionned conjugate gradient method. Return code:
  >0 upon completion, <0 no convergence
  0 memory problem 
 -1 ill conditioned matrix
 -2 max it reached  */
int csrPrecondGrad(pCsr A,double *x,double *b,double *er,int *ni,char tgv) {
  pCsr      L;
  double   *ap,*p,*q,*y;
  double    dp,nn,rmp,rm,rm2,alpha,beta,err;
  int       n,it,ier,nit;

  /* allocation */
  if ( !x || !b )  return(0);
  assert(A->nr == A->nc);
  n = A->nr;
  y = (double*)malloc(A->nr*sizeof(double));
  assert(y);

  /* compute R0 = y = b - A.x0 */
  nn = csrXY(x,x,A->nr);
  if ( nn < CS_EPSD2 ) {
    memcpy(y,b,A->nr*sizeof(double));
	}
	else {
    csrAxpy(A,x,b,y,-1.,1.);
  }
  rmp = csrXY(y,y,A->nr);
	if ( fabs(rmp) < CS_EPSD2 ) {
    free(y);
    return(1);
  }
  else if ( tgv ) {
    rmp /= CS_TGV2;
  }

  /* P_1 = P^-1.R_0 */
  p  = (double*)malloc(A->nr*sizeof(double));
  q  = (double*)malloc(A->nr*sizeof(double));
  ap = (double*)malloc(A->nr*sizeof(double));
  assert(p);
  assert(q);
  assert(ap);

  /* incomplete SSOR factorization */
  L = csrTr(A);
  csrSSOR(A,L,p,y);

  err = *er;
  err = err * err * rmp;
  nit = *ni;
  ier = 1;
  it  = 0;
  rm  = rmp;
  while ( (err < rm) && ++it <= nit ) {
    /* alpha_m = <P^-1.R_m-1,R_m-1> / <AP_m,P_m> */
    rm = csrXY(p,y,n);

		if ( fabs(rm) <= CS_EPSD2 )  break;
    dp = csrAxdotx(A,p,ap);
    if ( fabs(dp) <= CS_EPSD2 )  break;

    /* X_m = X_m-1 + alpha_m.P_m , R_m = R_m-1 - alpha_m AP_m */ 
    alpha = (rm / dp);
		csrlXmY(p,x,x,alpha,1.0,n);
		csrlXmY(ap,y,y,-alpha,1.0,n);

    /* beta_m = <P^-1.R_m,R_m> / <P^-1.R_m-1,R_m-1> */
    csrSSOR(A,L,q,y);
    rm2 = csrXY(q,y,n);
    if ( fabs(rm2) <= CS_EPSD2 )  break;

    /* P_m+1 = P^-1.R_m + beta_m P_m */
    beta = rm2 / rm;
    csrlXmY(p,q,p,beta,1.0,n);
    //printf("  GC1: it  %d err %E   rm %E\n",it,err,rm2); 
    rm = rm2;
  }
  if ( it > nit )   ier = -2;
  *er = sqrt(rm / rmp);
  *ni = it;
  free(p);
  free(q);
  free(y);
  free(ap);
  csrFree(L);
  return(ier);
}

/* Solve Ax=b using GMRES  Return code:
   >0 upon completion, <0 no convergence
   0 memory problem 
  -1 ill conditioned matrix
  -2 max it reached */
int csrGMRES(pCsr A,double *x,double *b,double *er,int *ni,int krylov,int prec) {
  double  *r,*V,h[CS_KRY1+1],y[CS_KRYL],g[CS_KRYL+1],c[CS_KRYL],s[CS_KRYL];
  double   beta,tol,dd,err;
  int      i,j,ier,restart,m,iadr,nit,it;

  if ( !x || !b )  return(0);
  assert(A->nr > 0);

  restart = krylov ? CS_MAX(krylov,CS_KRYL) : CS_KRYL;
  r = (double *)malloc(A->nr*sizeof(double));
  assert(r);

  /* compute R0 = b - A.x0 */
  memcpy(r,b,A->nr*sizeof(double));

  /* diagonal pre-cond */
  beta = 0.0;
  for (i=0; i<A->nr; i++) {
    if ( prec ) {
      for (j=A->row[i]; j<A->row[i+1]; j++) {
        if ( i == A->col[j] ) {
          assert(A->val[j]);
          r[i] /= A->val[j];
					break;
        }
      }
    }
    beta += r[i] * r[i];
  }
  if ( beta < CS_EPSD2 ) {
    fprintf(stdout,"  ## GMRES Pb: beta NULL %E\n",beta);
    free(r);
    return(0);
  }
  else if ( beta >= CS_TGV2 ) {
    beta /= CS_TGV2;
  }
  beta = sqrt(beta);
  dd   = 1.0 / beta;
  csrAxpy(A,x,r,NULL,-1.,1.);

  /* V1 = R0 / beta */
  V = (double *)malloc(A->nr*(restart+1)*sizeof(double));
	assert(V);
  for (i=0; i<A->nr; i++)  V[i] = r[i] * dd;

  it  = 0;
  nit = *ni;
  err = *er;
  tol = beta * err;
  ier = 1;
  while ( tol < beta && it < nit ) {
    m    = 0;
    iadr = 0;
    g[0] = beta;

    while ( (tol < fabs(g[m])) && (m < restart) ) {
      g[m+1] = 0;
      /* hi,m = (AVm, Vi) , i=1,..,m */
      csrAx(A,&V[m*A->nr],&V[(m+1)*A->nr]);

      if ( prec ) {
        for (i=0; i<A->nr; i++) {
          for (j=A->row[i]; j<A->row[i+1]; j++) {
            if ( i == A->col[j] ) {
              V[(m+1)*A->nr+i] /= A->val[j];
              break;
            }
          }
        }
      }

      for (i=m; i>=0; i--) {
        h[iadr+i] = 0.;
        for (j=0; j<A->nr; j++)
          h[iadr+i] += V[(m+1)*A->nr+j] * V[i*A->nr+j];
      }

      /* Vm+1 = AVm - sum(0,m, hi,m.Vi)  */    
      for (j=0; j<=m; j++) {
        for (i=0; i<A->nr; i++)
          V[(m+1)*A->nr+i] -= h[iadr+j] * V[j*A->nr+i];
      }

      /* hm+1,m = || Vm+1 || */
      h[iadr+m+1] = 0.0;
      for (i=0; i<A->nr; i++) 
        h[iadr+m+1] += V[(m+1)*A->nr+i] * V[(m+1)*A->nr+i]; 
      h[iadr+m+1] = sqrt(h[iadr+m+1]);
   
      if ( h[iadr+m+1] >  CS_EPSD2 ) {
        /* Vm+1 = Vm+1 / hm+1,m */
        dd = 1.0 / h[iadr+m+1];
        for (i=0; i<A->nr; i++)  V[(m+1)*A->nr+i] *= dd;
      }
      else
        fprintf(stdout,"  ## GMRES Pb: H NULL\n");

      /* update H */
      for (i=0; i<m; i++) {
        beta        = h[iadr+i];
        h[iadr+i]   = c[i]*h[iadr+i]   + s[i]*h[iadr+i+1];
        h[iadr+i+1] = c[i]*h[iadr+i+1] - s[i]*beta;
      }
      beta = sqrt(h[iadr+m]*h[iadr+m] + h[iadr+m+1]*h[iadr+m+1]);
      if ( beta < CS_EPSD2 ) {
        fprintf(stdout, "  ## GMRES Pb: beta NULL\n");
        free(r);  free(V);
        return(0);
      }
      dd   = 1.0 / beta;
      c[m] = h[iadr+m]   * dd;
      s[m] = h[iadr+m+1] * dd;

      beta = h[iadr+m];
      h[iadr+m]   = c[m]*h[iadr+m]   + s[m]*h[iadr+m+1];
      h[iadr+m+1] = c[m]*h[iadr+m+1] - s[m]*beta;

      /* update g */
      g[m+1] =  -1.0 * s[m]*g[m];
      g[m]   = c[m]*g[m];

      it++;
      m++;
      iadr = iadr + m + 1;
			if ( iadr > CS_KRY1 )  break;
    }
 
    /* compute y */
    iadr = iadr - m - 1;
    m--;
    for (j=m; j>=0; j--) {
      y[j] = g[j] / h[iadr+j];
      for (i=0; i<j; i++)
        g[i] -= h[iadr+i] * y[j];
      iadr = iadr -j -1;
    }

    /* Compute Un */  
    for (j=0; j<=m; j++) {
      for (i=0; i<A->nr; i++)
        x[i] += y[j] * V[j*A->nr+i];
    }    

    if ( it < nit && tol < fabs(g[m+1]) ) {
      /* Compute Un */      
      memcpy(r,b,A->nr*sizeof(double));
      csrAxpy(A,x,r,NULL,-1.,1.);

			beta = 0.0;
      for (i=0; i<A->nr; i++) {
        for (j=A->row[i]; j<A->row[i+1]; j++) {
          if ( prec && (i == A->col[j]) ) {
            assert(fabs(A->val[j]) > CS_EPSD2);
            r[i] /= A->val[j];
          }
        }
				beta += r[i] * r[i];
      }
      if ( beta < CS_EPSD2 ) {
        fprintf(stdout,"  ## GMRES Pb: Beta NULL\n");
        free(r);  free(V);
        return(0);
      }
      dd   = 1.0 / sqrt(beta);
      for (i=0; i<A->nr; i++)  V[i] = r[i] * dd;
    }
    beta = fabs(g[m+1]);
  }

  /* Memory free */
  free(r);  free(V);
  if ( it == nit )   ier = -2;
  *er = beta;
  *ni = it;
  return(ier);
}



