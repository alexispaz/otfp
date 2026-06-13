#include <stdio.h>
#include <stdlib.h>

#include "chapeau_internal.h"
#include "chapeau.h"

chapeau * chapeau_crop (chapeau * ch, double * rmin, double * rmax) {
  //Retrun in cho a portion (in z space) of ch
  int i,j;
  int min[2],max[2];
  int ni,nj,nk;
  int mi,mj,mk;
  chapeau * cho;

  // int * periodic[ch->dm],N[ch->dm];
  int * N;
  int * periodic;

  periodic=calloc(ch->dm,sizeof(int));
  N=calloc(ch->dm,sizeof(int));

  for (i=0;i<ch->dm;i++) {

    // if ( rmax[i] > ch->rmax[i] ) rmax[i]=ch->rmax[i];
    // if ( rmin[i] < ch->rmin[i] ) rmin[i]=ch->rmin[i];

    min[i]=(int)((rmin[i]-ch->rmin[i])*ch->idr[i]);
    max[i]=(int)((rmax[i]-ch->rmin[i])*ch->idr[i])+1;

    rmin[i]=min[i]*ch->dr[i]+ch->rmin[i];
    rmax[i]=max[i]*ch->dr[i]+ch->rmin[i];
    N[i]=max[i]-min[i]+1;

    // Lost periodicity... I guess
    periodic[i]=0;
  }


  cho=chapeau_alloc(ch->dm, rmin, rmax, N, periodic);

  for (i=0;i<ch->N[0];i++) {

    if (i<min[0]) continue;
    if (i>max[0]) continue;

    for (j=0;j<ch->N[1];j++) {
      if (j<min[1]) continue;
      if (j>max[1]) continue;

      // Identifico el punto de red y sus equivalentes en el nuevo espacio
      /*nj----------*
         |-\        |
         |  -\      |
         |    -\    |
         |  x   -\  |
         |        -\|
        nk---------ni */
      mi=j*ch->N[0]+i+1;
      mj=mi+ch->N[0]-1;
      mk=mi-1;

      ni=(j-min[1])*cho->N[0]+(i-min[0])+1;
      nj=ni+cho->N[0]-1;
      nk=ni-1;

      cho->b[nk] = ch->b[mk];
      cho->bfull[nk] = ch->bfull[mk];
      cho->hits[nk] = ch->hits[mk];
      cho->lam[nk] = ch->lam[mk];

      // Diagonal
      cho->A[cho->ku][nk] = ch->A[ch->ku][mk];  // This is A[nk][nk]
      cho->Afull[cho->ku][nk] = ch->Afull[ch->ku][mk];  // This is Afull[nk][nk]

      // Off diagonal
      if (i<max[0]) {
        cho->A[cho->ku+(ni-nk)][nk] = ch->A[ch->ku+(mi-mk)][mk]; // This is A[ni][nk]
        cho->A[cho->ku+(nk-ni)][ni] = ch->A[ch->ku+(mk-mi)][mi]; // This is A[nk][ni]
        cho->Afull[cho->ku+(ni-nk)][nk] = ch->Afull[ch->ku+(mi-mk)][mk]; // This is Afull[ni][nk]
        cho->Afull[cho->ku+(nk-ni)][ni] = ch->Afull[ch->ku+(mk-mi)][mi]; // This is Afull[nk][ni]
      }
      if (j<max[1]) {
        cho->A[cho->ku+(nj-nk)][nk] = ch->A[ch->ku+(mj-mk)][mk]; // This is A[nj][nk]
        cho->A[cho->ku+(nk-nj)][nj] = ch->A[ch->ku+(mk-mj)][mj]; // This is A[nk][nj]
        cho->Afull[cho->ku+(nj-nk)][nk] = ch->Afull[ch->ku+(mj-mk)][mk]; // This is Afull[nj][nk]
        cho->Afull[cho->ku+(nk-nj)][nj] = ch->Afull[ch->ku+(mk-mj)][mj]; // This is Afull[nk][nj]
      }
    }
  }

  // Border gross correction... FIXME
  // if (min[0]>0) {
  //   i=0;
  //   for (j=0;j<cho->N[1];j++) {
  //     ni=j*cho->N[0]+i;
  //     cho->A[cho->ku][ni]     = .5*cho->A[cho->ku][ni]     ;
  //     cho->Afull[cho->ku][ni] = .5*cho->Afull[cho->ku][ni] ;
  //     cho->b[ni]     = .5*cho->b[ni];
  //     cho->bfull[ni] = .5*cho->bfull[ni];
  //   }
  // }
  // if (max[0]<ch->N[0]) {
  //   i=cho->N[0];
  //   for (j=0;j<cho->N[1];j++) {
  //     ni=j*cho->N[0]+i;
  //     cho->A[cho->ku][ni]     = .5*cho->A[cho->ku][ni]     ;
  //     cho->Afull[cho->ku][ni] = .5*cho->Afull[cho->ku][ni] ;
  //     cho->b[ni]     = .5*cho->b[ni];
  //     cho->bfull[ni] = .5*cho->bfull[ni];
  //   }
  // }

  // This makes segfault
  // //Korner
  // j=max[1];
  // i=max[0];
  // mi=j*ch->N[0]+i+1;
  // mj=mi+ch->N[0]-1;
  // mk=mj+1;
  //
  // ni=(j-min[1])*cho->N[0]+(i-min[0])+1;
  // nj=ni+cho->N[0]-1;
  // nk=nj+1;
  //
  // cho->b[nk] = ch->b[mk];
  // cho->bfull[nk] = ch->bfull[mk];
  // cho->hits[nk] = ch->hits[mk];
  // cho->lam[nk] = ch->lam[mk];
  //
  // // Diagonal
  // cho->A[cho->ku][nk] = ch->A[ch->ku][mk];  // This is A[nk][nk]
  // cho->Afull[cho->ku][nk] = ch->Afull[ch->ku][mk];  // This is Afull[nk][nk]
  //
  // // Off diagonal
  // cho->A[cho->ku+(ni-nk)][nk] = ch->A[ch->ku+(mi-mk)][mk]; // This is A[ni][nk]
  // cho->A[cho->ku+(nk-ni)][ni] = ch->A[ch->ku+(mk-mi)][mi]; // This is A[nk][ni]
  // cho->Afull[cho->ku+(ni-nk)][nk] = ch->Afull[ch->ku+(mi-mk)][mk]; // This is Afull[ni][nk]
  // cho->Afull[cho->ku+(nk-ni)][ni] = ch->Afull[ch->ku+(mk-mi)][mi]; // This is Afull[nk][ni]
  // cho->A[cho->ku+(nj-nk)][nk] = ch->A[ch->ku+(mj-mk)][mk]; // This is A[nj][nk]
  // cho->A[cho->ku+(nk-nj)][nj] = ch->A[ch->ku+(mk-mj)][mj]; // This is A[nk][nj]
  // cho->Afull[cho->ku+(nj-nk)][nk] = ch->Afull[ch->ku+(mj-mk)][mk]; // This is Afull[nj][nk]
  // cho->Afull[cho->ku+(nk-nj)][nj] = ch->Afull[ch->ku+(mk-mj)][mj]; // This is Afull[nk][nj]

  return cho;
}

void chapeau_sum ( chapeau * ch1, chapeau * ch2 ) {
  //Overwrite ch1 with ch1+ch2
  int i,j;

  if (!chapeau_comparesize(ch1,ch2)) {
    fprintf(stderr,"OTFP: ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }

  //Really?
  //if (!chapeau_comparegrid(ch1,ch2)) {
  //  fprintf(stderr,"OTFP: ERROR: you can not sum chapeau objects with different domains\n");
  //  exit(-1);
  //}

  for (i=0;i<ch1->m;i++) {
    ch1->b[i] += ch2->b[i];
    ch1->bfull[i] += ch2->bfull[i];
    ch1->hits[i] = ch1->hits[i]+ch2->hits[i];
    for (j=0;j<ch1->ldad;j++) {
      ch1->A[j][i] += ch2->A[j][i];
      ch1->Afull[j][i] += ch2->Afull[j][i];
    }
  }

  //gsl_vector_add(ch1->lam ,ch2->lam); // irrelevant?

}

int chapeau_simetrize2D( chapeau * ch ) {
  int i,j,ii,jj,n,nn,m,mm;

  // Copy statistics across the line y=x
  // Por ahora solo funciona cuando N[0]=N[1] TODO: facil de corregir
  // No funciona con periodicidad (obvias razones, es imposible esta siemtria en ese caso)


  /* Recorro el triangulo inferior de la simetria

           +----------------+
           |              -/|
           |            -/  |
           |          -/    |
           |        -/      |
           |      -/        |
           |    -/          |
           |  -/     X      |
           |-/              |
           +----------------+

  */
  for (i=0;i<ch->N[0];i++) {
    for (j=0;j<i+1;j++) {

      // convierto la notacion
      n=j*ch->N[0]+i;

      // consigo su reflejo
      ii=j;
      jj=i;
      nn=jj*ch->N[0]+ii;

      // b y hits
      ch->b[n]              += ch->b[nn];
      ch->hits[n]            =(ch->hits[n]||ch->hits[nn]);

      // Diagonal central de A
      ch->A[ch->ku][n]      +=ch->A[ch->ku][nn];

      // Elementos fuera de la diagonal

      if(i!=ch->N[0]-1){

        /* 1st upper diagonal: Considerando el triangulo inverior de la simetrica
           esto es la interaccion horizontal con el nodo posterior, que en realidad
           en el triangulo superior es la interaccion con el nodo superior (Nth
           diagonal).

           +----------------+
           |  nn+1        -/|
           |   |        -/  |
           |   nn     -/    |
           |        -/      |
           |      -/        |
           |    -/          |
           |  -/   n-n+1    |
           |-/              |
           +----------------+

        */

        m =n+1;
        mm=nn+ch->N[0];

        ch->A[ch->ku+(n-m)][m]+=ch->A[ch->ku+(nn-mm)][mm];     // This is A[n][m]
        ch->A[ch->ku+(m-n)][n]+=ch->A[ch->ku+(mm-nn)][nn];     // This is A[m][n]
      }

      if(j!=i){

        /* Nst upper diagonal: Considerando el triangulo inverior de la simetrica
           esto es la interaccion horizontal con el nodo posterior, que en realidad
           en el triangulo superior es la interaccion con el nodo posterior (1st
           diagonal).

           +----------------+
           |              -/|
           |  nn-nn+1   -/  |
           |          -/    |
           |        -/      |
           |      -/   n+1  |
           |    -/      |   |
           |  -/        n   |
           |-/              |
           +----------------+

        */

        m =n+ch->N[0];
        mm=nn+1;

        ch->A[ch->ku+(n-m)][m]+=ch->A[ch->ku+(nn-mm)][mm];     // This is A[n][m]
        ch->A[ch->ku+(m-n)][n]+=ch->A[ch->ku+(mm-nn)][nn];     // This is A[m][n]
      }
      // Las interacciones hacia la izquierda y hacia abajo ya son contadas
      // en la interaccion anterior

    }
  }

  // Igualo al triangulo superior
  for (i=0;i<ch->N[0];i++) {
    for (j=i;j<ch->N[1];j++) {

      // convierto la notacion
      n=j*ch->N[0]+i;

      // consigo su reflejo
      ii=j;
      jj=i;
      nn=jj*ch->N[0]+ii;

      ch->b[n]         = ch->b[nn];
      ch->hits[n]      = (ch->hits[n]||ch->hits[nn]);
      ch->A[ch->ku][n] = ch->A[ch->ku][nn];

      if(i!=j){
        m = n+1;
        mm = nn+ch->N[0];
        ch->A[ch->ku+(n-m)][m] = ch->A[ch->ku+(nn-mm)][mm];     // This is A[n][m]
        ch->A[ch->ku+(m-n)][n] = ch->A[ch->ku+(mm-nn)][nn];     // This is A[m][n]
      }

      if(j!=ch->N[1]-1){
        m  = n+ch->N[0];
        mm = nn+1;
        ch->A[ch->ku+(n-m)][m] = ch->A[ch->ku+(nn-mm)][mm];     // This is A[n][m]
        ch->A[ch->ku+(m-n)][n] = ch->A[ch->ku+(mm-nn)][nn];     // This is A[m][n]
      }
    }
  }

  return 0;
}
