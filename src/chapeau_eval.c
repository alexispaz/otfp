#include <stdio.h>
#include <stdlib.h>

#include "chapeau_internal.h"
#include "chapeau.h"

#define SMALL 1.e-10

void chapeau_set_peaks ( chapeau * ch, char * filename ) {
  FILE * fp = fopen(filename,"r");
  double knots;
  int i=0;
  char ln[255];

  fprintf(stdout,"INFO) Read knots from %s\n",filename);

  for (i=0;i<ch->m;i++) {
    fgets(ln,255,fp);
    sscanf(ln,"%lf",&knots);
    fprintf(stdout,"INFO) %i,%.5f\n",i++,knots);
    ch->lam[i]=knots;
  }

  fflush(stdout);
}

int chapeau_caneval_f1 ( chapeau * ch, double z ) {
  int m;

  if (ch->dm!=1) {
    fprintf(stderr,"OTFP: ERROR: chapeau dimension does not match evaluation requested");
    exit(-1);
  }
  if ( z > ch->rmax[0] ) return 0;
  if ( z <= ch->rmin[0] ) return 0;

  m=(int) ch->idr[0]*(z-ch->rmin[0]);

  if (!ch->hits[m]) return 0;
  if (!ch->hits[m+1]) return 0;

  return 1;
}

int chapeau_caneval_f2 ( chapeau * ch, double z1, double z2 ) {
  int i,j,ni,nj,nk;
  double dx,dy;

  if (ch->dm!=2) {
    fprintf(stderr,"OTFP: ERROR: chapeau dimension does not match evaluation requested");
    exit(-1);
  }
  if ( z1 > ch->rmax[0] ) return 0;
  if ( z1 <= ch->rmin[0] ) return 0;
  if ( z2 > ch->rmax[1] ) return 0;
  if ( z2 <= ch->rmin[1] ) return 0;

  dx=(z1-ch->rmin[0])*ch->idr[0];
  dy=(z2-ch->rmin[1])*ch->idr[1];

  i=(int)(dx);
  j=(int)(dy);
  ni=j*ch->N[0]+i+1;
  nj=ni+ch->N[0]-1;

  dx=dx-i;
  dy=dy-j;
  if(dy<1-dx) {
    nk=ni-1;
  } else {
    nk=nj+1;
  }

  if (!ch->hits[ni]) return 0;
  if (!ch->hits[nj]) return 0;
  if (!ch->hits[nk]) return 0;

  return 1;
}

double chapeau_evalf_1simplex ( chapeau * ch, double z ) {
  int m;
  double dm,la,lb;

  dm=ch->idr[0]*(z-ch->rmin[0]);
  m=(int) dm;
  dm=dm-m;

  la=ch->lam[m];
  lb=ch->lam[m+1];
  return la+(lb-la)*dm;
}

double chapeau_evalf_2simplex ( chapeau * ch, double z1, double z2 ) {
  int i,j,ni,nj,nk;
  double f,dx,dy;

  dx=(z1-ch->rmin[0])*ch->idr[0];
  dy=(z2-ch->rmin[1])*ch->idr[1];

  i=(int)(dx);
  j=(int)(dy);
  ni=j*ch->N[0]+i+1;
  nj=ni+ch->N[0]-1;

  dx=dx-i;
  dy=dy-j;
  if(dy<1-dx) {
    nk=ni-1;

    f=(1+dx+dy)*ch->lam[nk];
    f+=dx*ch->lam[ni];
    f+=dy*ch->lam[nj];

  } else {
    nk=nj+1;

    f=(dx+dy-1)*ch->lam[nk];
    f+=(1-dy)*ch->lam[ni];
    f+=(1-dx)*ch->lam[nj];
  }

  return f;
}

double chapeau_evalfg_2simplex ( chapeau * ch, double z1, double z2, double * g ) {
  int i,j,ni,nj,nk;
  double f,dx,dy;

  dx=(z1-ch->rmin[0])*ch->idr[0];
  dy=(z2-ch->rmin[1])*ch->idr[1];

  if(z1==ch->rmax[0]) z1-=SMALL;
  if(z2==ch->rmax[1]) z2-=SMALL;

  i=(int)(dx);
  j=(int)(dy);
  ni=j*ch->N[0]+i+1;
  nj=ni+ch->N[0]-1;

  dx=dx-i;
  dy=dy-j;
  if(dy<1-dx) {
    nk=ni-1;

    f=(1-dx-dy)*ch->lam[nk];
    f+=dx*ch->lam[ni];
    f+=dy*ch->lam[nj];

    g[0]=-ch->lam[nk];
    g[0]+=ch->lam[ni];
    g[1]=-ch->lam[nk];
    g[1]+=ch->lam[nj];

  } else {
    nk=nj+1;

    f=(dx+dy-1)*ch->lam[nk];
    f+=(1-dy)*ch->lam[ni];
    f+=(1-dx)*ch->lam[nj];

    g[0]=ch->lam[nk];
    g[0]-=ch->lam[nj];
    g[1]=ch->lam[nk];
    g[1]-=ch->lam[ni];
  }

  g[0]=g[0]*ch->idr[0];
  g[1]=g[1]*ch->idr[1];

  return f;
}
