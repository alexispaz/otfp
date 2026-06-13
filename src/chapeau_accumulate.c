#include <math.h>

#include "chapeau_obj.h"
#include "chapeau.h"

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
