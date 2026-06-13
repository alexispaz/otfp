#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapeau_obj.h"
#include "chapeau.h"

int chapeau_init ( chapeau * ch, int dm, double * rmin, double * rmax, int * N, int * periodic) {
  int i;

  // IF NO DSS and (periodic[0]||periodic[1]) {
  //   fprintf(stderr,"Warning OTFP: Intel MKL Direct Sparse Solver is better for periodicty.\n");
  // }
  // IF DSS and !(periodic[0]||periodic[1]) {
  //   fprintf(stderr,"OTFP: No periodicity implies a banded matrix. Consider using lapack for better performance.\n");
  // }  
  
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


 int accumulate_1D( chapeau * ch, double bias ) { 
  int m;
  double boost;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_solve procedure)
  
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] > ch->rmax[0] ) return 0;
  if ( ch->r[0] <= ch->rmin[0] ) return 0;
 
  // Insted of solve Ax=b, I rather solve (dr[0]*dr[0]*A)x/dr[0]=(dr[0]*b).
  // Then, I have to remember multiply the solution by dr[0]. The factor idr[0]
  // in the terms of A and b are now factor 1. 
                       
  //ch->nsample++;

  //boost=1.;
  boost=exp(bias);

  /* Interpolation R->m
        m R
        | |
    o---o-x-o---o--*/
  m=(int)((ch->r[0]-ch->rmin[0])*ch->idr[0]);
  ch->hits[m]=1;
  ch->hits[m+1]=1; // Testing

  /* 1/dz  /| 
          / |n=m+1
         /  | 
    o---o-x-o---o--*/ 
  ch->b[m+1]+=ch->f[0]*boost;
  ch->A[ch->ku][m+1]+=boost;

  /*    |\   -1/dz
     n=m| \
        |  \
    o---o-x-o---o--*/ 
  ch->b[m]-=ch->f[0]*boost;
  ch->A[ch->ku][m]+=boost;

  /*    |\ /| 
       m| X |n=m+1
        |/ \| 
    o---o-x-o---o--*/ 
  ch->A[ch->ku+1][m]-=boost;
  ch->A[ch->ku-1][m+1]-=boost;


  if(ch->periodic[0]){
    if(m==ch->N[0]-2){

      // A hit in the rightmost knot => copy to leftmost knot
      ch->b[0]+=ch->f[0]*boost;
      ch->A[ch->ku][0]+=boost;
      ch->hits[0]=1;

      // Rightmost is actually not used at the moment of invert A (gives the
      // same equation of the leftmost).
      ch->hits[ch->N[0]-1]=0;
    }
   
  }

  return 0;
}
 
int accumulate_2D( chapeau * ch, double bias ) { 
  int i,j,ni,nj,nk,mi,mj,mk;
  double dx,dy,ratio,ratio2;
  double boost;
  
  // Following the paper E definition A and b should be
  //
  // A=1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz dphi_n(z)/dz
  // b=-1./nsteps \sum_{steps} \sum_{i} dphi_m(z)/dz F
  //
  // But here we compute -b*nsteps and A*nsteps and then A and b are scaled
  // by the proper factor when equation A\lambda=b is solved (see the
  // chapeau_solve procedure)
 
  // Early return to avoid interpolations out of the domain
  if ( ch->r[0] >= ch->rmax[0] ) return 0;
  if ( ch->r[1] >= ch->rmax[1] ) return 0;
  if ( ch->r[0] < ch->rmin[0] ) return 0;
  if ( ch->r[1] < ch->rmin[1] ) return 0;

  /* LA BASE
    
  Li, Z., Qiao, Z., & Tang, T. (2015). The Finite Element Method for 2D
  Problems. In Numerical Solutions of Partial differential equations- An
  introduction to finite difference and finite element methods (pp. 199–245).
  Url: http://www4.ncsu.edu/~zhilin/TEACHING/MA587/

  -\  |      -\  |      -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R5  |  -\      |  -\
      |    -\    |    -\    |    
  -\  |  R6  -\  | R4   -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\  R1  |  -\  R3  |  -\
      |    -\    |    -\    |    
  -\  |      -\  |  R2  -\  |    
    -\|        -\|        -\|    
  ----*----------*----------*----
      |-\        |-\        |-\
      |  -\      |  -\      |  -\
  
  */
             
  //boost=1.;
  boost=exp(bias);
            
  // The ratio allow me to transform the matrix A and vector b to make them
  // indpendent of the space. So, insted of solve Ax=b, I rather solve
  // (dr[0]*dr[0]*A)x/dr[0]=(dr[0]*b). Then, I have to remember multiply the
  // solution by dr[0]. The factor idr[0] in the terms of A and b are now
  // factor 1 and the terms with idr[1] are now the factor ratio. 
  ratio=ch->dr[0]*ch->idr[1]*boost;
  ratio2=ratio*ratio*boost;

  // Identifico el cuadrado, esto define dos nodos
  /*nj----------*
     |          |  
     |          |  
     |     x    |  
     |          |  
     |          |  
     *---------ni */
  dx=(ch->r[0]-ch->rmin[0])*ch->idr[0]; 
  dy=(ch->r[1]-ch->rmin[1])*ch->idr[1]; 

  i=(int)(dx);
  j=(int)(dy);
  ni=j*ch->N[0]+i+1;
  //nj=(j+1)*ch->N[0]+i;
  nj=ni+ch->N[0]-1;


  // The band storage scheme for a 
  // M = N = 6, KL = 2, KU = 1 example:
  //
  //   *    *    *    +    +    +      
  //   *    *    +    +    +    + 
  //   *   a01  a12  a23  a34  a45 (uperdiagonal)
  //  a00  a11  a22  a33  a44  a55 (diagonal)
  //  a10  a21  a32  a43  a54   * 
  //  a20  a31  a42  a53   *    * 
  //
  // storate[(kl+ku+(i-j))][j]=a[i][j]; 
  //
  // Note that columns index are the same
 
           
  // ahora el triangulo, define el tercer nodo
  dx=dx-i; 
  dy=dy-j;
  if(dy<1-dx) {

    //nk=j*ch->N[0]+i;
    nk=ni-1;

    /*nj----------* 
       |-\        | 
       |  -\      | 
       |    -\    | 
       |  x   -\  | 
       |        -\| 
      nk---------ni */
    
    ch->b[ni]         += ch->f[0]*boost;
    ch->A[ch->ku][ni] += boost;            // This is A[ni][ni]

    ch->b[nj]         += ratio*ch->f[1];
    ch->A[ch->ku][nj] += ratio2;         // This is A[nj][nj]
     
    ch->b[nk]        -= ch->f[0]*boost;
    ch->b[nk]        -= ratio*ch->f[1];
    ch->A[ch->ku][nk]+= boost;             // This is A[nk][nk]
    ch->A[ch->ku][nk]+= ratio2;          // This is A[nk][nk]
    
    // TODO: Aca podria evitar 2 operaciones si eligiera el
    // triangulo superior de A (aprovechando la simetria)
    ch->A[ch->ku+(ni-nk)][nk]-= boost;     // This is A[ni][nk]
    ch->A[ch->ku+(nk-ni)][ni]-= boost;     // This is A[nk][ni]
    ch->A[ch->ku+(nj-nk)][nk]-= ratio2;  // This is A[nj][nk]
    ch->A[ch->ku+(nk-nj)][nj]-= ratio2;  // This is A[nk][nj]

    ch->hits[ni]=1;
    ch->hits[nj]=1;
    ch->hits[nk]=1;
 
    if(ch->periodic[0]){
      if(i==ch->N[0]-2){

        // A hit in the rightmost knot => copy to leftmost knot
        mi=ni-ch->N[0]+1;
        ch->b[mi]         += ch->f[0]*boost;
        ch->A[ch->ku][mi] += boost;
        ch->hits[mi]=1;

        // Rightmost is actually not used at the moment of invert A (gives
        // the same equation of the leftmost).
        ch->hits[ni]=0;

        // Interaction across boundary
        ch->A[ch->ku+(mi-nk)][nk]-= boost;     // This is A[ni][nk]
        ch->A[ch->ku+(nk-mi)][mi]-= boost;     // This is A[nk][ni]
      }
    }  
    if(ch->periodic[1]){
      if(j==ch->N[1]-2){

        // A hit in the uppermost knot => copy to lowermost knot
        mj=i;
        ch->b[mj]         += ratio*ch->f[1];
        ch->A[ch->ku][mj] += ratio2;
        ch->hits[mj]=1;

        // Lowermost is actually not used at the moment of invert A (gives
        // the same equation of the uppermost).
        ch->hits[nj]=0; 

        // Interaction across boundary
        ch->A[ch->ku+(mj-nk)][nk] -= ratio2;  // This is A[nj][nk] 
        ch->A[ch->ku+(nk-mj)][mj] -= ratio2;  // This is A[nk][nj] 
      }
    }
  } else {
    //nk=(j+1)*ch->N[0]+ch->ku+1;
    nk=nj+1;
 
    /*nj---------nk
       |-\        |  
       |  -\   x  |  
       |    -\    |  
       |      -\  |  
       |        -\|  
       *---------ni */

    ch->b[ni]         -= ratio*ch->f[1];
    ch->A[ch->ku][ni] += ratio2;
                     
    ch->b[nj]         -= ch->f[0]*boost;
    ch->A[ch->ku][nj] += boost;
                     
    ch->b[nk]         += ch->f[0]*boost;
    ch->b[nk]         += ratio*ch->f[1];
    ch->A[ch->ku][nk] += boost;
    ch->A[ch->ku][nk] += ratio2;

    // TODO: Aca podria evitar 2 operaciones si eligiera el
    // triangulo superior de A (aprovechando la simetria)
    ch->A[ch->ku+(nk-nj)][nj] -= boost;      // This is A[nk][nj]
    ch->A[ch->ku+(nj-nk)][nk] -= boost;      // This is A[nj][nk]
    ch->A[ch->ku+(nk-ni)][ni] -= ratio2;   // This is A[nk][ni]
    ch->A[ch->ku+(ni-nk)][nk] -= ratio2;   // This is A[ni][nk]
 
    ch->hits[ni]=1;
    ch->hits[nj]=1;
    ch->hits[nk]=1;

    if(ch->periodic[0]){
      if(i==ch->N[0]-2){

        // A hit in the rightmost knot => copy to leftmost knot
        mi=ni-ch->N[0]+1;
        ch->b[mi]         -= ratio*ch->f[1];
        ch->A[ch->ku][mi] += ratio2;
        ch->hits[mi]=1;

        mk=nk-ch->N[0]+1;
        ch->b[mk]         += ch->f[0]*boost;
        ch->b[mk]         += ratio*ch->f[1];
        ch->A[ch->ku][mk] += boost;
        ch->A[ch->ku][mk] += ratio2;
        ch->hits[mk]=1;
        
        // Terminos cruzados que sobreviven (los otros estan en la otra replica)
        ch->A[ch->ku+(mk-mi)][mi] -= ratio2;   // This is A[nk][ni]
        ch->A[ch->ku+(mi-mk)][mk] -= ratio2;   // This is A[ni][nk]

        // Rightmost is actually not used at the moment of invert A (gives
        // the same equation of the leftmost).
        ch->hits[ni]=0;
        ch->hits[nk]=0;

        // Interaction across boundary
        ch->A[ch->ku+(mk-nj)][nj] -= boost;    // This is A[nk][nj]
        ch->A[ch->ku+(nj-mk)][mk] -= boost;    // This is A[nj][nk]

        // Esta esta fuera de alcance, porque la fila nk no se usa
        // ch->A[ch->ku+(nk-mi)][mi] -= ratio2;   // This is A[nk][ni]
        // ch->A[ch->ku+(mi-nk)][mk] -= ratio2;   // This is A[ni][nk]

      }

    }  
    if(ch->periodic[1]){
      if(j==ch->N[1]-2){
        // hit in the uppermost knot => copy to lowermost knot
        mj=i;
        ch->b[mj]         -= ch->f[0]*boost;
        ch->A[ch->ku][mj] += boost;
        ch->hits[mj]=1;
              
        mk=mj+1;
        ch->b[mk]         += ch->f[0]*boost;
        ch->b[mk]         += ratio*ch->f[1];
        ch->A[ch->ku][mk] += boost;
        ch->A[ch->ku][mk] += ratio2;
        ch->hits[mk]=1;
        
        // terminos cruzados que sobreviven (los otros estan en la otra replica)
        ch->A[ch->ku+(mk-mj)][mj] -= boost;      // This is A[nk][nj]
        ch->A[ch->ku+(mj-mk)][mk] -= boost;      // This is A[nj][nk]
 
        // Uppermost is actually not used at the moment of invert A (gives
        // the same equation of the lowermost).
        ch->hits[nj]=0; 
        ch->hits[nk]=0; 

        // Interaction across boundary
        ch->A[ch->ku+(mk-ni)][ni] -= ratio2;   // This is A[nk][ni]
        ch->A[ch->ku+(ni-mk)][mk] -= ratio2;   // This is A[ni][nk]
 
      }
    }

    if(ch->periodic[0] && ch->periodic[1]){
      if(i==ch->N[0]-2 && j==ch->N[1]-2){

        // A hit the upper right korner => copy to the lower left korner  
        mk=0;
        ch->b[mk]         += ch->f[0]*boost;
        ch->b[mk]         += ratio*ch->f[1];
        ch->A[ch->ku][mk] += boost;
        ch->A[ch->ku][mk] += ratio2;
        ch->hits[mk]=1;   

        // Esta esta fuera de alcance, porque la fila ni no se usa
        // ch->A[ch->ku+(mk-ni)][ni] -= ratio2;   // This is A[nk][ni]
        // ch->A[ch->ku+(ni-mk)][mk] -= ratio2;   // This is A[ni][nk]
        // ch->A[ch->ku+(mk-nj)][nj] -= boost;    // This is A[nk][nj]
        // ch->A[ch->ku+(nj-mk)][mk] -= boost;    // This is A[nj][nk]
      }
    }
  }
  return 0;
}
 
