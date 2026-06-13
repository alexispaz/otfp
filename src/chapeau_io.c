#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chapeau_obj.h"
#include "chapeau.h"

void chapeau_setupoutput ( chapeau * ch, char * outfile, char * restartfile, int outputFreq, int outputLevel ) {

  if (!ch) exit(-1);

  // output
  ch->ofp=fopen(outfile,"w");
  char * flag="CxCx";
  fwrite(flag,sizeof(char),4,ch->ofp);
  fwrite(&outputLevel,sizeof(int),1,ch->ofp);
  fwrite(&ch->m,sizeof(int),1,ch->ofp);
  ch->outputFreq=outputFreq;
  ch->outputLevel=outputLevel;

  // restart
  strcpy(ch->filename,restartfile);

}

void chapeau_output ( chapeau * ch, int timestep ) {
  int outputlevel=ch->outputLevel;
  int i;

  if ( !ch || !ch->ofp || timestep%ch->outputFreq ) return;

  fwrite(&timestep,sizeof(int),1,ch->ofp);

  if (outputlevel & 1) { // 0th bit = output knots as y
    for (i=0;i<ch->m;i++) {
      fwrite(&(ch->lam[i]),sizeof(double),1,ch->ofp);
    }
  }

  if (outputlevel & 2) {
    for (i=0;i<ch->m;i++) {
      fwrite(&(ch->hits[i]),sizeof(int),1,ch->ofp);
    }
  }

  fflush(ch->ofp);
}

void chapeau_savestate ( chapeau * ch, char * filename ) {
  FILE *ofs;
  int i;

  ofs=fopen(filename,"w");

  fwrite(ch, sizeof(*ch), 1, ofs);

  //TODO
  // fwrite(ch->periodic.. );
  // fwrite(ch->norm.. );
  // fwrite(ch->ku.. );

  fwrite(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
  fwrite(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
  fwrite(ch->N,sizeof(*(ch->N)),ch->dm,ofs);

  // getchar();
  fwrite(ch->lam,sizeof(*(ch->lam)),ch->m,ofs);
  // getchar();
  fwrite(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);

  for (i=0;i<ch->ldad;i++) {
    fwrite(ch->A[i],sizeof(*(ch->A[i])),ch->m,ofs);
  }
  fwrite(ch->b,sizeof(*(ch->b)),ch->m,ofs);

  for (i=0;i<ch->ldad;i++) {
    fwrite(ch->Afull[i],sizeof(*(ch->Afull[i])),ch->m,ofs);
  }
  fwrite(ch->bfull,sizeof(*(ch->bfull)),ch->m,ofs);

  fclose(ofs);

}

void chapeau_loadstate ( chapeau * ch, char * filename ) {
    chapeau * chaux;
    int i;
    FILE * ofs;

    ofs=fopen(filename,"r");

    if (!ch) {
      fprintf(stderr,"CFACV) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }

    // With chaux we prevent override of addresses (hits, A, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    chaux = malloc(sizeof(chapeau));
    fread(chaux, sizeof(*chaux), 1,ofs);
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"CFACV) ERROR, the chapeau object from the file is diferent\n");
    }

    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m    = chaux->m;
    ch->ldad    = chaux->ldad;

    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;

    // Read directly
    fread(ch->rmin,sizeof(*(ch->rmin)),ch->dm,ofs);
    fread(ch->rmax,sizeof(*(ch->rmax)),ch->dm,ofs);
    fread(ch->N,sizeof(*(ch->N)),ch->dm,ofs);

    for (i=0;i<ch->dm;i++) {
      ch->dr[i]  =(ch->rmax[i]-ch->rmin[i])/(ch->N[i]-1);
      ch->idr[i] =1./ch->dr[i];
    }

    fread(ch->lam,sizeof(*(ch->lam)),ch->m,ofs);
    fread(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->A[i],sizeof(*(ch->A[i])),ch->m,ofs);
    }
    fread(ch->b,sizeof(*(ch->b)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->Afull[i],sizeof(*(ch->Afull[i])),ch->m,ofs);
    }
    fread(ch->bfull,sizeof(*(ch->bfull)),ch->m,ofs);

    fclose(ofs);
    free(chaux);
}

chapeau * chapeau_allocloadstate ( char * filename ) {
    int i;
    FILE * ofs;
    chapeau * ch;
    double * rmin;
    double * rmax;
    int * N;
    int * periodic;

    fprintf(stdout,"CFACV) Allocating chapeau object from restart file\n");

    ofs=fopen(filename,"r");

    ch = malloc(sizeof(chapeau));

    // Read overall information of chapeau object
    fread(ch, sizeof(*ch), 1,ofs);

    // Alloc minimal info to run initializacion of the chapeau
    periodic=calloc(ch->dm,sizeof(int));
    rmin=malloc(ch->dm*sizeof(double));
    rmax=malloc(ch->dm*sizeof(double));
    N=malloc(ch->dm*sizeof(int));

    // Read the minimal informacion
    // TODO: Read periodic vector
    fread(rmin,sizeof(*rmin),ch->dm,ofs);
    fread(rmax,sizeof(*rmax),ch->dm,ofs);
    fread(N,sizeof(*N),ch->dm,ofs);

    // Initialize object
    chapeau_init(ch,ch->dm,rmin,rmax,N,periodic);

    // Read the rest of the information
    fread(ch->lam,sizeof(*(ch->lam)),ch->m,ofs);
    fread(ch->hits,sizeof(*(ch->hits)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->A[i],sizeof(*(ch->A[i])),ch->m,ofs);
    }
    fread(ch->b,sizeof(*(ch->b)),ch->m,ofs);

    for (i=0;i<ch->ldad;i++) {
      fread(ch->Afull[i],sizeof(*(ch->Afull[i])),ch->m,ofs);
    }
    fread(ch->bfull,sizeof(*(ch->bfull)),ch->m,ofs);

    fclose(ofs);

    return ch;
}

void chapeau_loadlambda ( chapeau * ch, char * filename ) {
    int i;
    chapeau * chaux;

    if (!ch) {
      fprintf(stderr,"OTFP) ERROR in load chapeau file because holding object was not allocated\n");
      exit(-1);
    }

    chaux=chapeau_allocloadstate(filename);

    // With chaux we prevent override of addresses (hits, A, etc) when
    // reading integers (N,M,et), that is only what we want... better way to do
    // this?
    if (!chapeau_comparesize(ch,chaux)) {
      fprintf(stdout,"OTFP) ERROR, the chapeau object from the file is diferent\n");
    }

    // This variables are alredy set during allocation, might be better to
    // compare them instead
    ch->dm      = chaux->dm;
    ch->m    = chaux->m;
    ch->ldad    = chaux->ldad;

    for (i=0;i<ch->dm;i++) {
      ch->rmin[i] = chaux->rmin[i];
      ch->rmax[i] = chaux->rmax[i];
      ch->N[i] = chaux->N[i];
      ch->dr[i] = (ch->rmax[i]-ch->rmin[i])/(ch->m-1);
      ch->idr[i] = 1./ch->dr[i];
    }

    // This variables might be not set
    ch->nupdate     = chaux->nupdate    ;
    strcpy(ch->filename,chaux->filename);
    ch->outputFreq  = chaux->outputFreq ;
    ch->outputLevel = chaux->outputLevel;

    for (i=0;i<ch->m;i++) {
      ch->lam[i]=chaux->lam[i];
    }

    chapeau_free(chaux);
}

char * chapeau_serialize ( chapeau * ch ) {
  // return a str containing the serialized chapeau with partial statistics. If
  // this is the central replica, or non replica scheme is used, this serialize
  // the full statistic.
  int i, j;
  int size;
  char* buf;

  // size = (3*ch->m-1)*13+ch->m+1;
  size  = ch->m*13;           // space for ch->b
  size += ch->m*ch->ldad*13;; // space for ch->A
  size += ch->m;;             // space for ch->hits
  size += 1;;                 // null character?
  buf  = malloc(size*sizeof(char));

  // Seems that if I do this:
  //  return (char*)&ch
  //  is not portable, since the way that the cast is made is undefined
  //  (depends of the machine). So it is needed a proper serialization

  size=0;
  for (i=0;i<ch->m;i++) {

    // Writing the partial b
    sprintf(buf+size, "%13.5e", ch->b[i]); size+=13;

    // Writing hits
    sprintf(buf+size, "%1i", ch->hits[i]); size+=1;
  }

  // TODO take advantage of the simetric property to reduce half of the numbers
  // Writing the partial A (banded)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      sprintf(buf+size, "%13.5e",ch->A[i][j]); size+=13;
    }
  }

  return buf;
}

void chapeau_addserialized ( chapeau *ch, char * str ) {
  // str contains the serialized chapeau that comes from the partial sampling
  // of other replica partial statistic information. This is added to the
  // partial sampling of the principal replica computing the sum.
  int i,j,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];

  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';

  size1  = ch->m*13;           // space for ch->b
  size1 += ch->m*ch->ldad*13;; // space for ch->A
  size1 += ch->m;;             // space for ch->hits
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"OTFP: ERROR: you can not sum chapeaus with different sizes %d %d\n",size1,size2);
    exit(-1);
  }

  // Warning, sscanf, strtod, strtoi, all start in the first nonblank
  // position, then there is not a proper way to read the serialized object
  // using sscanf.
  //
  // TODO, strtod and strtoi more efficient that sscanf?

  size1=0;
  for (i=0;i<ch->m;i++) {

    // Reading the partial b
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"OTFP: Error 1112 on read %s\n",word1);}
    ch->b[i]+=aux;

    // Reading hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"OTFP: Error 1111 on read %s\n",word2);}
    ch->hits[i]=(ch->hits[i]||iaux);
  }

  // Reading the partial A (simetric and tridiagonal)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      memcpy(word1, &str[size1], 13 ); size1+=13;
      err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"OTFP: Error 1112 on read %s\n",word1);}
      ch->A[i][j]+=aux;
    }
  }

  //TODO, give error if rmin and rmax are different

}

void chapeau_setserialized ( chapeau *ch, char * str ) {
  // str contains the full statistics information that comes from all the
  // replicas contributions. Therefore, this should be stored in Afull and
  // bfull and the partial A and b should be reset to cero.
  int i,j,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];

  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';

  size1  = ch->m*13;           // space for ch->b
  size1 += ch->m*ch->ldad*13;; // space for ch->A
  size1 += ch->m;;             // space for ch->hits
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"OTFP: ERROR: you can not set chapeaus with different sizes %d %d\n",size1,size2);
    exit(-1);
  }

  // Warning, sscanf, strtod, strtoi, all start in the first nonblank
  // position, then there is not a proper way to read the serialized object
  // using sscanf.
  //
  // TODO, strtod and strtoi more efficient that sscanf?

  size1=0;
  for (i=0;i<ch->m;i++) {

    // Read b
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"OTFP: Error 2111 on read %s\n",word1);}
    ch->bfull[i]=aux;

    // Read hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"OTFP: Error 2112 on read %s\n",word1);}
    ch->hits[i]=iaux;
  }

  // Reading the full A (simetric and tridiagonal)
  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) {
      memcpy(word1, &str[size1], 13 ); size1+=13;
      err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"OTFP: Error 1112 on read %s\n",word1);}
      ch->Afull[i][j]=aux;
    }
  }

  // Reset the variables to hold the partial statistic information
  for (j=0;j<ch->m;j++) ch->b[j]=0.;

  for (i=0;i<ch->ldad;i++) {
    for (j=0;j<ch->m;j++) ch->A[i][j]=0.;
  }

  //TODO, give error if rmin and rmax are different
}
