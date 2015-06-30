#include "chapeau.h"

void chapeau_sum ( chapeau * ch1, chapeau * ch2 ) {
  //Overwrite ch1 with ch1+ch2
  int i,j,k;

  if (sizeof(ch1)!=sizeof(ch2)) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }
  if (ch1->rmin!=ch2->rmin || ch1->rmax!=ch2->rmax || ch1->dr!=ch2->dr ) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different sizes\n");
    exit(-1);
  }

  for (i=0;i<ch1->N;i++) {
    if ( ch1->mask[i]!=ch2->mask[i] ) {
      fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeau objects with different masks\n");
      exit(-1);
    }
  }    

  gsl_vector_add(ch1->b   ,ch2->b);
  gsl_matrix_add(ch1->A   ,ch2->A);
  gsl_vector_add(ch1->bfull,ch2->bfull);
  gsl_matrix_add(ch1->Afull,ch2->Afull);

  for (i=0;i<ch1->m;i++) {
    ch1->hits[i] = (ch1->hits[i]||ch2->hits[i]);
  }

  //// I think that the next is irrelevant since I am saving the state before
  //// the chapeau_output and therefore before the chapeau_update, but I let this
  //// here just for the case.
  //for (i=0;i<ch1->N;i++) {
  //  for (j=0;j<3;j++) {
  //    for (k=0;k<ch1->m;k++) {
  //      ch1->s[i][j][k] = ch1->s[i][j][k] + ch2->s[i][j][k];
  //    }
  //  }
  //}
  //gsl_vector_add(ch1->lam ,ch2->lam); //also irrelevant

}

chapeau * chapeau_alloc ( int m, double rmin, double rmax, int npart ) {
  chapeau * ch = (chapeau*)malloc(sizeof(chapeau));
  int d,i;

  if (npart<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without particles\n");
    exit(-1);
  }         

  if (m<1) {
    fprintf(stderr,"CFACV/C) ERROR: chapeau objects without bins\n");
    exit(-1);
  }         

  ch->m=m;
  ch->N=npart;
  ch->rmin=rmin;
  ch->rmax=rmax;
  ch->dr=(rmax-rmin)/(m-1);
  ch->idr=1.0/ch->dr;
  ch->lam=gsl_vector_calloc(m);
  ch->nsample=0;
  ch->hits=(int*)calloc(m,sizeof(int));
  
  //ch->s=(double***) calloc(npart,sizeof(double**));
  //for (i=0;i<npart;i++) {
  //  ch->s[i]=(double**)calloc(3,sizeof(double*));
  //  for (d=0;d<3;d++) {
  //    ch->s[i][d]=(double*)calloc(m,sizeof(double));
  //  }
  //}

  ch->mask=(int*)calloc(npart,sizeof(int));
  
  ch->b=gsl_vector_calloc(m);
  ch->A=gsl_matrix_calloc(m,m);
  ch->bfull=gsl_vector_calloc(m);
  ch->Afull=gsl_matrix_calloc(m,m);

  fprintf(stdout,"CFACG4/DEBUG) allocated chapeau with %i knots on [%.5f,%.5f] dr %g \n",m,rmin,rmax,ch->dr);
  fflush(stdout);

  ch->alpha=0;

  return ch;
}

void chapeau_setPeaks ( chapeau * ch, double * peaks ) {
  if (ch) {
    if (peaks) {
      int i;
      for (i=0;i<ch->m;i++) gsl_vector_set(ch->lam,i,peaks[i]);
    }
  }
}


void chapeau_setupoutput ( chapeau * ch, char * filename, int outputFreq, int outputLevel ) {
  char foo[255];
    
  if (ch) {
    char * flag="CxCx";
    ch->ofp=fopen(filename,"w");

    //header
    fwrite(flag,sizeof(char),4,ch->ofp);
    fwrite(&outputLevel,sizeof(int),1,ch->ofp);
    fwrite(&ch->m,sizeof(int),1,ch->ofp);

    ch->outputFreq=outputFreq;
    ch->outputLevel=outputLevel;

    strcpy(foo,filename);
    strcat(foo,".restart");
    ch->ofs=fopen(foo,"w");
  }
}

void chapeau_output ( chapeau * ch, int timestep ) {
  if (ch&&ch->ofp&&!(timestep%ch->outputFreq)) {
    int outputlevel=ch->outputLevel;
    int i;
  
    fwrite(&timestep,sizeof(int),1,ch->ofp);
    if (outputlevel & 1) { // 0th bit = output knots as y
      for (i=0;i<ch->m;i++) {
	fwrite(gsl_vector_ptr(ch->lam,i),sizeof(double),1,ch->ofp);
	//fprintf(stderr,"### %i %g\n",i,*gsl_vector_ptr(ch->lam,i));
      }
    }
    if (outputlevel & 2) {
      for (i=0;i<ch->m;i++) {
	fwrite(&(ch->hits[i]),sizeof(int),1,ch->ofp);
      }
    }
    fflush(ch->ofp);

    chapeau_savestate(ch,timestep);
  }
}

void chapeau_savestate ( chapeau * ch, int timestep ) {
  if (ch&&ch->ofs&&!(timestep%ch->outputFreq)) {
    int N=ch->N;
    int m=ch->m;

    rewind(ch->ofs);
    fwrite(&timestep,sizeof(int),1,ch->ofs);

    fwrite(ch, sizeof(*ch), 1, ch->ofs);
    fwrite(ch->hits,sizeof(*(ch->hits)),m,ch->ofs);
    fwrite(ch->mask,sizeof(*(ch->mask)),N,ch->ofs);
    //fwrite(ch->s,sizeof(***(ch->s)),m*N*3,ch->ofs);

    gsl_matrix_fwrite(ch->ofs,ch->A);
    gsl_vector_fwrite(ch->ofs,ch->b);
    gsl_matrix_fwrite(ch->ofs,ch->Afull);
    gsl_vector_fwrite(ch->ofs,ch->bfull);
    gsl_vector_fwrite(ch->ofs,ch->lam);
    fflush(ch->ofs);
  }
}

chapeau * chapeau_allocloadstate ( char * filename ) {
    int i,m,N,d;
    FILE * ofs=fopen(filename,"r");
    chapeau * ch;
   
    fread(&i,sizeof(int),1,ofs);
    fprintf(stdout,"## Recovering chapeau state of step %i of some previous simulation\n",i);

    ch = (chapeau*)malloc(sizeof(chapeau));
  
    fread(ch, sizeof(*ch), 1, ofs);
    N=ch->N;
    m=ch->m;
   
    ch->hits=(int*)malloc(m*sizeof(int));
    fread(ch->hits,sizeof(*(ch->hits)),m,ofs);

    ch->mask=(int*)malloc(N*sizeof(int));
    fread(ch->mask,sizeof(*(ch->mask)),N,ofs);

    //ch->s=(double***) malloc(N*sizeof(double**));
    //for (i=0;i<N;i++) {
    //  ch->s[i]=(double**)malloc(3*sizeof(double*));
    //  for (d=0;d<3;d++) {
    //    ch->s[i][d]=(double*)malloc(m*sizeof(double));
    //    fread(ch->s[i][d],sizeof(*(ch->s[i][d])),m,ofs);
    //  }
    //}

    ch->A=gsl_matrix_calloc(m,m);
    gsl_matrix_fread(ofs,ch->A);
    ch->b=gsl_vector_calloc(m);
    gsl_vector_fread(ofs,ch->b);

    ch->Afull=gsl_matrix_calloc(m,m);
    gsl_matrix_fread(ofs,ch->Afull);
    ch->bfull=gsl_vector_calloc(m);
    gsl_vector_fread(ofs,ch->bfull);

    ch->lam=gsl_vector_calloc(m);
    gsl_vector_fread(ofs,ch->lam);

    fread(&ch,sizeof(*ch),1,ofs);
    fclose(ofs);
 
    return ch;
}
    
void chapeau_update_peaks ( chapeau * ch ) {
  int i,j,J,I;
  int s;
  double ninv;
  double lo,lb,alpha;

  gsl_matrix * Abar;
  gsl_vector * bbar;
  gsl_vector * lambar;
  gsl_permutation * p;
  int nred;

  //DB//for (i=0;i<ch->m;i++) fprintf(stderr,"AAA %i %i\n",i,ch->hits[i]); exit(1);

  //// parameters that are allowed to evolve lie between indices for which
  //// ch->hits[] is non-zero so extract the proper subspace
  nred=0;
  for (i=0;i<ch->m;i++) if (ch->hits[i]) nred++;

  // If the nred is small, singular matrix can occur
  //fprintf(stderr,"Warning: No chapeau update: only %i non cero elements\n",nred);
  if (nred<10) return;

  Abar=gsl_matrix_alloc(nred,nred);
  bbar=gsl_vector_alloc(nred);
  lambar=gsl_vector_alloc(nred);
  p=gsl_permutation_alloc(nred);
  

  // Add the new and old statistic to the reduced matrix
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;
    lb=gsl_vector_get(ch->b,i)+gsl_vector_get(ch->bfull,i);
    gsl_vector_set(bbar,I,-lb); //TODO: Trace back the origin of the minus. See fes1D procedure. 
    J=0;
    for (j=0;j<ch->m;j++) {
      if (!ch->hits[j]) continue;
      lb=gsl_matrix_get(ch->A,i,j)+gsl_matrix_get(ch->Afull,i,j);
      gsl_matrix_set(Abar,I,J,lb);
      J++;
    }
    I++;
  } 

  //XXX: What is alpha?
  alpha=0.0;//1.0-exp(-1.0e4/nsamples);

  gsl_linalg_LU_decomp(Abar,p,&s);
  gsl_linalg_LU_solve(Abar,p,bbar,lambar);
  
  // update the vector of coefficients
  I=0;
  for (i=0;i<ch->m;i++) {
    if (!ch->hits[i]) continue;

    lb=gsl_vector_get(lambar,I);
    if (lb!=lb) {fprintf(stderr,"PARANOIA) Tripped at 3333. Too many chapeau additions?\n");exit(-1);}
    gsl_vector_set(ch->lam,i,lb);
    
    //XXX: What is alpha?
    //lo=gsl_vector_get(ch->lam,i);
    //gsl_vector_set(ch->lam,i,alpha*lo+(1-alpha)*lb);

    I++;
  } 
   
  gsl_matrix_free(Abar);
  gsl_vector_free(bbar);
  gsl_vector_free(lambar);
  gsl_permutation_free(p);
}


// Process for replica exechange
    
double chapeau_evalf ( chapeau * ch, double z ) {
  // Evaluate the free energy of a given CV vector from the linear
  // interpolation of the lambda vector on a chapeau object.
  // For now the interpolation is 1D.
  int m;
  double dm,la,lb;

  // Early return to avoid interpolations beyond the boundaries
  if ( z > ch->rmax ) return 0.;
  if ( z <= ch->rmin ) return 0.;
   
  // o                    z   
  // |                o   |   
  // |       o        |   |   o
  // |       |        |   |   |
  // |       |        |   |   |
  // min  min+dr ...  |   v   |
  // +-------+---//---+-------+---
  // 0       1   ...  m      m+1
  // /---------dm---------/

  dm=ch->idr*(z-ch->rmin);

  m=(int) dm;
  dm=dm-m;
  la=gsl_vector_get(ch->lam,m);
  lb=gsl_vector_get(ch->lam,m+1);
  return la+(lb-la)*dm; 

}  
    
char * chapeau_serialize ( chapeau * ch ) {
  // return a str containing the serialized chapeau with partial statistics. If
  // this is the central replica, or non replica scheme is used, this serialize
  // the full statistic.
  int i;
  int size;
  char* buf;

  size = (3*ch->m-1)*13+ch->m+1;
  buf  = (char*) malloc(size*sizeof(char));

  // Seems that if I do this:
  //  return (char*)&ch
  //  is not portable, since the way that the cast is made is undefined
  //  (depends of the machine). So it is needed a proper serialization
  
  size=0;
  for (i=0;i<ch->m;i++) {

    // Writing the partial b
    sprintf(buf+size, "%13.5e", gsl_vector_get(ch->b,i)); size+=13;

    // Writing hits
    sprintf(buf+size, "%1i", ch->hits[i]); size+=1;
  }

  // Writing the partial A (simetric and tridiagonal)
  sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,0,0)); size+=13;
  for (i=1;i<ch->m;i++) {
    sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,i,i)); size+=13;
    sprintf(buf+size, "%13.5e",gsl_matrix_get(ch->A,i,i-1)); size+=13; 
  }

  return buf; 
}     

void chapeau_addserialized ( chapeau *ch, char * str ) {
  // str contains the serialized chapeau that comes from the partial sampling
  // of other replica partial statistic information. This is added to the
  // partial sampling ot the principal replica computing the sum.
  int i,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';

  size1 = (3*ch->m-1)*13+ch->m;
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"CFACV/C) ERROR: you can not sum chapeaus with different sizes %d %d\n",size1,size2);
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
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
    aux+=gsl_vector_get(ch->b,i);
    gsl_vector_set(ch->b,i,aux);

    // Reading hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 1111 on read %s\n",word2);}
    ch->hits[i]=(ch->hits[i]||iaux);
  }
 
  // Reading the partial A (simetric and tridiagonal)
  memcpy(word1, &str[size1], 13 ); size1+=13;
  err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1112 on read %s\n",word1);}
  aux+=gsl_matrix_get(ch->A,0,0);
  gsl_matrix_set(ch->A,0,0,aux);

  for (i=1;i<ch->m;i++) {

    //diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1113 on read %s\n",word1);}
    aux+=gsl_matrix_get(ch->A,i,i);
    gsl_matrix_set(ch->A,i,i,aux);              

    //offdiagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 1114 on read %s\n",word1);}
    aux+=gsl_matrix_get(ch->A,i,i-1);
    gsl_matrix_set(ch->A,i,i-1,aux);
    gsl_matrix_set(ch->A,i-1,i,aux);
  }
   

  //TODO, give error if rmin and rmax are different
  
}

void chapeau_setserialized ( chapeau *ch, char * str ) {
  // str contains the full statistics information that comes from all the
  // replicas contributions. Therefore, this should be stored in Afull and
  // bfull and the partial A and b should be reset to cero.
  int i,err;
  int size1;
  int size2;
  double aux;
  int iaux;
  char word1[14],word2[2];
  
  // Add null terminators
  word1[13]='\0';
  word2[1]='\0';
  
  size1 = (3*ch->m-1)*13+ch->m;
  size2 = strlen(str);
  if (size1!=size2) {
    fprintf(stderr,"CFACV/C) ERROR: you can not set chapeaus with different sizes %d %d\n",size1,size2);
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
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2111 on read %s\n",word1);}
    gsl_vector_set(ch->bfull,i,aux);

    // Read hits
    memcpy(word2, &str[size1], 1 ); size1+=1;
    err=sscanf(word2,"%1i",&iaux); if(!err) {fprintf(stderr,"CFACV/C) Error 2112 on read %s\n",word1);}
    ch->hits[i]=iaux;
  }
     
  // Read A (simetric and tridiagonal)
  memcpy(word1, &str[size1], 13 ); size1+=13; 
  err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2113 on read %s\n",word1);}
  gsl_matrix_set(ch->Afull,0,0,aux);

  for (i=1;i<ch->m;i++) {

    //diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2114 on read %s\n",word1);}
    gsl_matrix_set(ch->Afull,i,i,aux);

    //off diagonal
    memcpy(word1, &str[size1], 13 ); size1+=13;
    err=sscanf(word1,"%13le",&aux); if(!err) {fprintf(stderr,"CFACV/C) Error 2115 on read %s\n",word1);}
    gsl_matrix_set(ch->Afull,i,i-1,aux);
    gsl_matrix_set(ch->Afull,i-1,i,aux);
  }
  
  // Reset the variables to hold the partial statistic information
  gsl_vector_set_zero(ch->b);
  gsl_matrix_set_zero(ch->A);

  //TODO, give error if rmin and rmax are different
}
 
