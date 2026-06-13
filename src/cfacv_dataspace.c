#include <stdio.h>
#include <stdlib.h>

#include "cfacv_internal.h"
#include "cfacv.h"
#include "chapeau_internal.h"
#include "cvs.h"
#include "measurements.h"

DataSpace * NewDataSpace ( int Nctr, int Ncv, int Nrstr, long int seed ) {
  DataSpace * ds = calloc(1,sizeof(DataSpace));
  int i;
  ds->N=Nctr;
  ds->K=Nrstr;
  ds->iN=0;
  ds->iK=0;

  ds->ncv=Ncv;   ds->icv=0;

  // This will not work with numbers less that 16 bits!
  //// Random seed. I deatached this from the data space a do it global
  //// for all the file (see for example HarmonicCart_cutoff function)
  //Xi=malloc(3*sizeof(unsigned short));
  //Xi[0]=(seed & 0xffff00000000) >> 32;
  //Xi[1]=(seed & 0x0000ffff0000) >> 16;
  //Xi[2]= 0x330e;
  srand48(seed);
  fprintf(stderr,"CFACV) first drand48 will be %.5f\n",drand48());
  fprintf(stderr,"CFACV) first gasdev  will be %.5f and %.5f\n",gasdev(),gasdev());
  srand48(seed);

  //pfout=fopen("/dev/urandom","r");
  //fread(&kkk,sizeof(long int),1,pfout);
  //srand48(kkk);
  //oldptr=seed48(&Xi[0]);

  ds->R=malloc(Nctr*sizeof(double*));
  for (i=0;i<Nctr;i++) ds->R[i]=calloc(3,sizeof(double));

  ds->ac = calloc(Nctr,sizeof(atomCenterStruct*));

  ds->cv=calloc(Ncv,sizeof(cv*));
  ds->restr=calloc(Nrstr,sizeof(restraint*));

  ds->doAnalyticalCalc=0;
  ds->evolveAnalyticalParameters=0;
  ds->beginaccum=0;

  ds->ch=NULL;

  sprintf(ds->filename, "checkpoint.rstr");
  ds->restrsavefreq = 100000;

  return ds;
}

void DataSpace_free ( DataSpace * ds ) {
  int i;

  if (!ds) return;
  if (ds->R) {
    for (i=0;i<ds->N;i++) free(ds->R[i]);
    free(ds->R);
  }
  if (ds->ac) {
    for (i=0;i<ds->N;i++) atomCenterStruct_free(ds->ac[i]);
    free(ds->ac);
  }
  if (ds->cv) {
    for (i=0;i<ds->ncv;i++) cv_free(ds->cv[i]);
    free(ds->cv);
  }
  if (ds->restr) {
    for (i=0;i<ds->K;i++) restraint_free(ds->restr[i]);
    free(ds->restr);
  }
  if (ds->ch) {
    for (i=0;i<ds->ch_num;i++) chapeau_free(ds->ch[i]);
    free(ds->ch);
  }
  free(ds);
}

chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i ) {
  return ds->ch[i];
}

int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double*  min, int * nKnots,
    double * max, int * periodic, int beginaccum, int beginsolve, char * outfile,
    int outfreq, int outlevel, int nupdate) {
  // Initialize chapeaus objects in DataSpace. Since Data space might not have
  // any chapeau allocated, the number of chapeaus is not indicated in the
  // dataspace constructor.
  int i;
  char buf[80], buf2[80];

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->doAnalyticalCalc=1;

  // This is to let the system equilibrate
  ds->evolveAnalyticalParameters=1;
  ds->beginaccum=beginaccum;
  ds->beginsolve=beginsolve;
  if (beginaccum < 0) ds->evolveAnalyticalParameters=0;

  ds->ch_num=numrep;
  ds->ch=malloc(ds->ch_num*sizeof(chapeau*));

  for (i=0;i<ds->ch_num;i++) {
    ds->ch[i]=chapeau_alloc(dm,min,max,nKnots,periodic);
    sprintf(buf, "%s_ch%d.bsp", outfile,i);
    sprintf(buf2, "%s.ch%d", outfile,i);
    chapeau_setupoutput(ds->ch[i],buf,buf2,outfreq,outlevel);
    chapeau_set_update_frequency(ds->ch[i],nupdate);
  }

  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  return ds->N;
}

double * DataSpace_centerPos ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->N) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return ds->R[i];
}

int DataSpace_checkdata ( DataSpace * ds ) {
  int i,j;

  if (!ds) return 1;

  for (i=0;i<ds->N;i++) {
    for (j=0;j<3;j++) {
      if (ds->R[i][j]!=ds->R[i][j]) {
        return 1;
      }
    }
  }

  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->N;i++) {
      fprintf(stdout,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->R[i][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
  return 0;
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  int i,j;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->iN++;
  if (ds->iN>ds->N) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->iN-1;
  ds->ac[i]=New_atomCenterStruct(n);
  for (j=0;j<n;j++) {
    ds->ac[i]->ind[j]=ind[j];
    ds->ac[i]->m[j]=m[j];
    ds->ac[i]->M+=m[j];
  }

  return i;
}

cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf,double boundk,
          char * outfile, int outputFreq )  {
  int i;
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->icv++;
  if (ds->icv>ds->ncv) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->icv-1;
  ds->cv[i]=cv_init(typ,nind,ind,zmin,zmax,boundf,boundk,outfile,outputFreq);

  return ds->cv[i];
}

restraint * DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr,
			 double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq )
{
  int i;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->iK++;
  if (ds->iK>ds->K) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->iK-1;
  ds->restr[i]=New_restraint(k,z,nCV,cvc,rftypstr,zmin,zmax,boundf,boundk,outfile, outputFreq);

  //restart
  if (outputFreq>0) {
    if (outputFreq < ds->restrsavefreq) {ds->restrsavefreq = outputFreq;}
  }

  return ds->restr[i];
}
