#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MKL_LAPACKE
  #include <mkl_lapacke.h>
#else
  #include <lapacke.h>
#endif

#include "chapeau_internal.h"
#include "chapeau.h"

void chapeau_solve ( chapeau * ch ) {
  int i,j,jj,ii,s,k;
  double * Abar;
  double * bbar;
  int * pivot;
  int nred;
  int ldad;
  char aux;

  nred=0;
  for (i=0;i<ch->m;i++)  if (ch->hits[i]) nred++;

  if (nred<ch->ku+1) {
   fprintf(stderr,"Warning: No chapeau update, nred (%d) < nro diagonals (%d)\n",nred,ch->ku);
   return;
  }

  fprintf(stdout,"OTFP) Inverting matrix A (%dx%d)\n",nred,nred);
  fflush(stdout);

  ldad=ch->ldad+ch->ku;
  bbar=malloc(nred*sizeof(double));
  Abar=calloc(nred*ldad,sizeof(double));
  pivot=malloc(nred*sizeof(int));

  jj=0;
  for (j=0;j<ch->m;j++) {
    if (!ch->hits[j]) continue;
    bbar[jj]=-(ch->b[j]+ch->bfull[j]);

    ii=jj;
    for (i=j;i<=j+ch->ku;i++) {
      if (i>=ch->m) break;
      if (!ch->hits[i]) continue;

      Abar[2*ch->ku+ii-jj+ldad*jj]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      ii++;
    }

    ii=jj-1;
    for (k=1;k<=ch->ku;k++) {
      i=j-k;
      if (i<0) break;
      if (!ch->hits[i]) continue;

      Abar[2*ch->ku+ii-jj+ldad*jj]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      ii--;
    }
    jj++;
  }

  if (nred!=jj) {fprintf(stderr,"Bad matrix nred size: %d != %d\n",nred,jj);exit(-1);}

  jj=1;
  aux='N';

  s = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, nred, nred, ch->ku, ch->ku, Abar, ldad, pivot);

  if (s!=0) {
    fprintf(stdout,"Matrix nearly singular with flag: %d\n",s);

    free(Abar);
    free(bbar);
    free(pivot);

    s=s-1;
    ii=0;
    for (i=0;i<ch->m;i++) {
      if (!ch->hits[i]) continue;
      if (ii==s) {
        fprintf(stdout,"Removing row %d\n",ii);
        ch->hits[i]=0;
        break;
      }
      ii++;
    }
    fflush(stdout);

    chapeau_solve(ch);

    return;
  }

  s = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, aux, nred, ch->ku, ch->ku, jj, Abar, ldad, pivot, bbar, nred);

  ii=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;
    ch->lam[i]=bbar[ii]*ch->dr[0];
    bbar[ii]=0;
    ii++;
  }

  free(Abar);
  free(bbar);
  free(pivot);
}

void chapeau_solve_secure ( chapeau * ch ) {
  int i,j,jj,ii,k,s;
  double lb;
  double * Abar;
  double * Ebar;
  double * Dbar;
  double * nullbar = NULL;
  double * U;
  double * VT;
  double * bbar;
  double * work;
  int nred;
  char aux;

  nred=0;
  for (i=0;i<ch->m;i++)  if (ch->hits[i]) nred++;

  if (nred<ch->ku+1) {
   fprintf(stderr,"Warning: No chapeau update: 0 columns (%d) > diagonals (%d)\n",nred,ch->ku+1);
   return;
  }

  bbar=malloc(nred*sizeof(double));
  work=malloc(4*nred*sizeof(double));
  Abar=calloc(nred*ch->ldad,sizeof(double));
  Dbar=calloc(nred,sizeof(double));
  Ebar=calloc(nred-1,sizeof(double));
  U=calloc(nred*nred,sizeof(double));
  VT=calloc(nred*nred,sizeof(double));

  jj=0;
  for (j=0;j<ch->m;j++) {
    if (!ch->hits[j]) continue;
    bbar[jj]=-(ch->b[j]+ch->bfull[j]);

    ii=jj;
    for (i=j;i<=j+ch->ku;i++) {
      if (i>=ch->m) break;
      if (!ch->hits[i]) continue;

      Abar[ch->ku+ii-jj+ch->ldad*jj]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      ii++;
    }

    ii=jj-1;
    for (k=1;k<=ch->ku;k++) {
      i=j-k;
      if (i<0) break;
      if (!ch->hits[i]) continue;

      Abar[ch->ku+ii-jj+ch->ldad*jj]=ch->A[ch->ku+i-j][j]+ch->Afull[ch->ku+i-j][j];
      ii--;
    }
    jj++;
  }

  if (nred!=jj) {fprintf(stderr,"Bad matrix nred size: %d != %d\n",nred,jj);exit(-1);}

  ii=0;
  jj=1;
  aux='B';

  s = LAPACKE_dgbbrd(LAPACK_COL_MAJOR, aux, nred, nred, ii, ch->ku, ch->ku, Abar, ch->ldad, Dbar, Ebar,
          U, nred, VT, nred, nullbar, jj);

  if (s!=0) {fprintf(stderr,"dgbbrd: %d argument had illegal value\n",s);exit(-1);}

  aux='U';
  s = LAPACKE_dbdsqr(LAPACK_COL_MAJOR, aux,nred, nred, nred, ii, Dbar, Ebar, VT, nred, U, nred, nullbar, jj);
  if (s!=0) {fprintf(stderr,"dbdsqr: error with flag %d\n",s);exit(-1);}

  ii=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;
    lb=0;

    for (j=0;j<nred;j++) {
      if (Dbar[j]<1.e-10) {
        lb+=0.;
      } else {
        for (k=0;k<nred;k++)  lb+=VT[j+nred*ii]*U[k+nred*j]*bbar[k]/Dbar[j];
      }
    }

    ch->lam[i]=lb*ch->dr[0];
    ii++;
  }

  free(Abar);
  free(Dbar);
  free(Ebar);
  free(VT);
  free(U);
  free(bbar);
  free(work);
}
