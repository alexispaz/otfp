#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chapeau_internal.h"
#include "chapeau.h"

int chapeau_init ( chapeau * ch, int dm, double * rmin, double * rmax, int * N, int * periodic) {
  int i;

  //TODO: check size of argument vectors and write the errors
  
  ch->ofp=NULL;
  ch->nupdate     = 0;
  ch->outputFreq  = 0;
  ch->outputLevel = 0;

  // Initialize filename bytes
  memset(ch->filename, 0, 255*sizeof(char));

  // Allocating the dimension
  ch->dm=dm;
  ch->periodic=calloc(dm,sizeof(int));
  ch->rmin=calloc(dm,sizeof(double));
  ch->rmax=calloc(dm,sizeof(double));
  ch->dr=calloc(dm,sizeof(double));
  ch->idr=calloc(dm,sizeof(double));
  ch->N=calloc(dm,sizeof(int));
  ch->r=calloc(dm,sizeof(double));
  ch->f=calloc(dm,sizeof(double));

  if (ch->dm==1) {
    ch->accumulate=accumulate_1D;
  } else {
    ch->accumulate=accumulate_2D;
  }

  ch->m=1;
  for (i=0;i<dm;i++) {
    ch->N[i]=N[i];   
    ch->rmin[i]=rmin[i];
    ch->rmax[i]=rmax[i];
    ch->dr[i]  =(rmax[i]-rmin[i])/(ch->N[i]-1);
    ch->idr[i] =1./ch->dr[i];
    ch->m=ch->m*N[i];
  }
          
  // Allocating the vectors
  ch->lam=calloc(ch->m,sizeof(double));
  ch->hits=calloc(ch->m,sizeof(int));
  ch->b=calloc(ch->m,sizeof(double));
  ch->bfull=calloc(ch->m,sizeof(double));
  ch->norm=1.;
                 
  // chapeau 1D => 1 uperdiagonals
  // chapeau 2D => ch->N[0] uperdiagonals
  if (ch->dm==1) {
    ch->ku=1;
  } else{
    // TODO: En realidad si entiendo bien, la mayoria de estas diagonales van a tener ceros.
    // salvo la primera y la ch->N[0], osea que aca quizas se podría dejar
    // en 2 no N[0] si es que se pudiera empaquetar así.
    ch->ku=ch->N[0];
  }  
  ch->ldad=2*ch->ku+1;

  // The stiffness matrix:
  // An m-by-n band matrix with kl subdiagonals and ku
  // superdiagonals may be stored compactly in a
  // two-dimensional array with kl+ku+1 rows and n columns.
  // Columns of the matrix are stored in corresponding columns
  // of the array, and diagonals of the matrix are stored in
  // rows of the array.
  // The band storage scheme for a m = n = 6, kl = 2, ku = 2 example:
  //
  //   *    *   a02  a13  a24  a35 (second uperdiagonal)
  //   *   a01  a12  a23  a34  a45 (first uperdiagonal)
  //  a00  a11  a22  a33  a44  a55 (diagonal)
  //  a10  a21  a32  a43  a54   *  (first lowerdiagonal)
  //  a20  a31  a42  a53   *    *  (second lowerdiagonal) 
  //
  // storage[ku+(i-j)][j]=a[i][j]; 
  //
  // TODO: This storage can be reduced if I choose to store
  // the upper or lower triangle since the matrix is
  // symmetric.
              
  ch->A=calloc(ch->ldad,sizeof(double*));
  ch->Afull=calloc(ch->ldad,sizeof(double*));
  for (i=0;i<ch->ldad;i++) {
    ch->A[i]=calloc(ch->m,sizeof(double));
    ch->Afull[i]=calloc(ch->m,sizeof(double));
  }

  return 0;
}
 
chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int * periodic) {
  chapeau * ch;

  // Allocating the object
  ch=malloc(sizeof(chapeau));
  chapeau_init(ch, dm, rmin, rmax, N, periodic);

  return ch;
}

void chapeau_free ( chapeau * ch ) {
  int i;

  if (!ch) return;
  if (ch->ofp) fclose(ch->ofp);
  free(ch->periodic);
  free(ch->rmin);
  free(ch->rmax);
  free(ch->dr);
  free(ch->idr);
  free(ch->N);
  free(ch->r);
  free(ch->f); 

  free(ch->lam);
  free(ch->hits);
  free(ch->b);
  free(ch->bfull);
  for (i = 0; i < ch->ldad; i++) { 
      free(ch->A[i]);
      free(ch->Afull[i]);
  }
  free(ch->A);
  free(ch->Afull);
  
  free(ch);
  ch=NULL;
}
             
int chapeau_comparesize ( chapeau * ch1,  chapeau * ch2) {
  if (ch1->dm       != ch2->dm      ) return 0;
  if (ch1->m        != ch2->m       ) return 0;
  if (ch1->ldad     != ch2->ldad    ) return 0;
  return 1;
}
 
int chapeau_comparegrid ( chapeau * ch1,  chapeau * ch2) {
  int i;
  //if (ch1->nupdate     != ch2->nupdate    ) return 0;
  //if (ch1->outputFreq  != ch2->outputFreq ) return 0;
  //if (ch1->outputLevel != ch2->outputLevel) return 0;
  //if (strcmp(ch1->filename,ch2->filename)) return 0;
  if (ch1->dm          != ch2->dm         ) return 0;
  for (i=0;i<ch2->dm;i++) {
    if (ch1->rmin[i]  != ch2->rmin[i] ) return 0;
    if (ch1->rmax[i]  != ch2->rmax[i] ) return 0;
    if (ch1->dr[i]    != ch2->dr[i]   ) return 0;
    if (ch1->idr[i]   != ch2->idr[i]  ) return 0;
    if (ch1->N[i]  != ch2->N[i] ) return 0;
  }
  return 1;
}

void chapeau_set_update_frequency ( chapeau * ch, int nupdate ) {
  ch->nupdate=nupdate;
}

void chapeau_set_sample ( chapeau * ch, int dim, double r, double f ) {
  ch->r[dim]=r;
  ch->f[dim]=f;
}

int chapeau_accumulate_sample ( chapeau * ch, double bias ) {
  return ch->accumulate(ch,bias);
}

int chapeau_should_solve ( chapeau * ch, int timestep, int beginsolve ) {
  return timestep>beginsolve && !(timestep % ch->nupdate);
}

int chapeau_should_output ( chapeau * ch, int timestep ) {
  return !(timestep % ch->outputFreq);
}

void chapeau_savestate_default ( chapeau * ch ) {
  chapeau_savestate(ch,ch->filename);
}
