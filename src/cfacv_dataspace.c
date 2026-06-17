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
  ds->center_capacity=Nctr;
  ds->restraint_capacity=Nrstr;
  ds->center_count=0;
  ds->restraint_count=0;

  ds->cv_capacity=Ncv;   ds->cv_count=0;

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

  ds->center_data=malloc(Nctr*sizeof(double*));
  for (i=0;i<Nctr;i++) ds->center_data[i]=calloc(3,sizeof(double));

  ds->centers = calloc(Nctr,sizeof(atomCenterStruct*));

  ds->cvs=calloc(Ncv,sizeof(cv*));
  ds->restraints=calloc(Nrstr,sizeof(restraint*));

  ds->doAnalyticalCalc=0;
  ds->evolveAnalyticalParameters=0;
  ds->beginaccum=0;

  ds->chapeaus=NULL;

  sprintf(ds->filename, "checkpoint.rstr");
  ds->restraint_save_frequency = 100000;

  return ds;
}

void DataSpace_free ( DataSpace * ds ) {
  int i;

  if (!ds) return;
  if (ds->center_data) {
    for (i=0;i<ds->center_capacity;i++) free(ds->center_data[i]);
    free(ds->center_data);
  }
  if (ds->centers) {
    for (i=0;i<ds->center_capacity;i++) atomCenterStruct_free(ds->centers[i]);
    free(ds->centers);
  }
  if (ds->cvs) {
    for (i=0;i<ds->cv_capacity;i++) cv_free(ds->cvs[i]);
    free(ds->cvs);
  }
  if (ds->restraints) {
    for (i=0;i<ds->restraint_capacity;i++) restraint_free(ds->restraints[i]);
    free(ds->restraints);
  }
  if (ds->chapeaus) {
    for (i=0;i<ds->chapeau_count;i++) chapeau_free(ds->chapeaus[i]);
    free(ds->chapeaus);
  }
  free(ds);
}

chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i ) {
  return ds->chapeaus[i];
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

  ds->chapeau_count=numrep;
  ds->chapeaus=malloc(ds->chapeau_count*sizeof(chapeau*));

  for (i=0;i<ds->chapeau_count;i++) {
    ds->chapeaus[i]=chapeau_alloc(dm,min,max,nKnots,periodic);
    sprintf(buf, "%s_ch%d.bsp", outfile,i);
    sprintf(buf2, "%s.ch%d", outfile,i);
    chapeau_setupoutput(ds->chapeaus[i],buf,buf2,outfreq,outlevel);
    chapeau_set_update_frequency(ds->chapeaus[i],nupdate);
  }

  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  return ds->center_capacity;
}

double * DataSpace_centerPos ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->center_capacity) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return ds->center_data[i];
}

int DataSpace_checkdata ( DataSpace * ds ) {
  int i,j;

  if (!ds) return 1;

  for (i=0;i<ds->center_capacity;i++) {
    for (j=0;j<3;j++) {
      if (ds->center_data[i][j]!=ds->center_data[i][j]) {
        return 1;
      }
    }
  }

  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->center_capacity;i++) {
      fprintf(stdout,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->center_data[i][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
  return 0;
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  int i,j;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->center_count++;
  if (ds->center_count>ds->center_capacity) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->center_count-1;
  ds->centers[i]=New_atomCenterStruct(n);
  for (j=0;j<n;j++) {
    ds->centers[i]->ind[j]=ind[j];
    ds->centers[i]->m[j]=m[j];
    ds->centers[i]->M+=m[j];
  }

  return i;
}

cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf,double boundk,
          char * outfile, int outputFreq )  {
  int i;
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->cv_count++;
  if (ds->cv_count>ds->cv_capacity) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->cv_count-1;
  ds->cvs[i]=cv_init(typ,nind,ind,zmin,zmax,boundf,boundk,outfile,outputFreq);

  return ds->cvs[i];
}

restraint * DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr,
			 double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq )
{
  int i;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->restraint_count++;
  if (ds->restraint_count>ds->restraint_capacity) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->restraint_count-1;
  ds->restraints[i]=New_restraint(k,z,nCV,cvc,rftypstr,zmin,zmax,boundf,boundk,outfile, outputFreq);

  //restart
  if (outputFreq>0) {
    if (outputFreq < ds->restraint_save_frequency) {ds->restraint_save_frequency = outputFreq;}
  }

  return ds->restraints[i];
}
