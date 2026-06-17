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
  ds->ctr_total=Nctr;
  ds->rstr_total=Nrstr;
  ds->ctr_count=0;
  ds->rstr_count=0;

  ds->cv_total=Ncv;   ds->cv_count=0;

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

  ds->ctr_positions=malloc(Nctr*sizeof(double*));
  ds->ctr_forces=malloc(Nctr*sizeof(double*));
  for (i=0;i<Nctr;i++) {
    ds->ctr_positions[i]=calloc(3,sizeof(double));
    ds->ctr_forces[i]=calloc(3,sizeof(double));
  }

  ds->ctrs = calloc(Nctr,sizeof(atomCenterStruct*));

  ds->cvs=calloc(Ncv,sizeof(cv*));
  ds->rstrs=calloc(Nrstr,sizeof(restraint*));

  ds->fes_enabled=0;
  ds->fes_evolve_parameters=0;
  ds->fes_begin_accum=0;

  ds->fes=NULL;

  sprintf(ds->filename, "checkpoint.rstr");
  ds->rstr_save_frequency = 100000;

  return ds;
}

void DataSpace_free ( DataSpace * ds ) {
  int i;

  if (!ds) return;
  if (ds->ctr_positions) {
    for (i=0;i<ds->ctr_total;i++) free(ds->ctr_positions[i]);
    free(ds->ctr_positions);
  }
  if (ds->ctr_forces) {
    for (i=0;i<ds->ctr_total;i++) free(ds->ctr_forces[i]);
    free(ds->ctr_forces);
  }
  if (ds->ctrs) {
    for (i=0;i<ds->ctr_total;i++) atomCenterStruct_free(ds->ctrs[i]);
    free(ds->ctrs);
  }
  if (ds->cvs) {
    for (i=0;i<ds->cv_total;i++) cv_free(ds->cvs[i]);
    free(ds->cvs);
  }
  if (ds->rstrs) {
    for (i=0;i<ds->rstr_total;i++) restraint_free(ds->rstrs[i]);
    free(ds->rstrs);
  }
  if (ds->fes) {
    for (i=0;i<ds->fes_total;i++) chapeau_free(ds->fes[i]);
    free(ds->fes);
  }
  free(ds);
}

chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i ) {
  return ds->fes[i];
}

int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double*  min, int * nKnots,
    double * max, int * periodic, int fes_begin_accum, int fes_begin_solve, char * outfile,
    int outfreq, int outlevel, int nupdate) {
  // Initialize FES state. The current implementation stores each FES state in
  // a chapeau object because chapeau is the active basis/kernel.
  int i;
  char buf[80], buf2[80];

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->fes_enabled=1;

  // This is to let the system equilibrate
  ds->fes_evolve_parameters=1;
  ds->fes_begin_accum=fes_begin_accum;
  ds->fes_begin_solve=fes_begin_solve;
  if (fes_begin_accum < 0) ds->fes_evolve_parameters=0;

  ds->fes_total=numrep;
  ds->fes=malloc(ds->fes_total*sizeof(chapeau*));

  for (i=0;i<ds->fes_total;i++) {
    ds->fes[i]=chapeau_alloc(dm,min,max,nKnots,periodic);
    sprintf(buf, "%s_ch%d.bsp", outfile,i);
    sprintf(buf2, "%s.ch%d", outfile,i);
    chapeau_setupoutput(ds->fes[i],buf,buf2,outfreq,outlevel);
    chapeau_set_update_frequency(ds->fes[i],nupdate);
  }

  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  return ds->ctr_total;
}

double * DataSpace_centerPosition ( DataSpace * ds, int i ) {
  return DataSpace_ctr_position(ds,i);
}

double * DataSpace_centerForce ( DataSpace * ds, int i ) {
  return DataSpace_ctr_force(ds,i);
}

double * DataSpace_ctr_position ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return ds->ctr_positions[i];
}

double * DataSpace_ctr_force ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return ds->ctr_forces[i];
}

void DataSpace_clear_ctr_forces ( DataSpace * ds ) {
  int i,d;

  for (i=0;i<ds->ctr_total;i++) {
    for (d=0;d<3;d++) ds->ctr_forces[i][d]=0.0;
  }
}

void DataSpace_add_ctr_force ( DataSpace * ds, int ctr, int dim, double force ) {
  ds->ctr_forces[ctr][dim]+=force;
}

int DataSpace_checkdata ( DataSpace * ds ) {
  int i,j;

  if (!ds) return 1;

  for (i=0;i<ds->ctr_total;i++) {
    for (j=0;j<3;j++) {
      if (ds->ctr_positions[i][j]!=ds->ctr_positions[i][j]) {
        return 1;
      }
      if (ds->ctr_forces[i][j]!=ds->ctr_forces[i][j]) {
        return 1;
      }
    }
  }

  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->ctr_total;i++) {
      fprintf(stdout,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->ctr_positions[i][j]);
      fprintf(stdout,"| ");
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->ctr_forces[i][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
  return 0;
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  int i,j;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->ctr_count++;
  if (ds->ctr_count>ds->ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->ctr_count-1;
  ds->ctrs[i]=New_atomCenterStruct(n);
  for (j=0;j<n;j++) {
    ds->ctrs[i]->ind[j]=ind[j];
    ds->ctrs[i]->m[j]=m[j];
    ds->ctrs[i]->M+=m[j];
  }

  return i;
}

cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf,double boundk,
          char * outfile, int outputFreq )  {
  int i;
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->cv_count++;
  if (ds->cv_count>ds->cv_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->cv_count-1;
  ds->cvs[i]=cv_init(typ,nind,ind,zmin,zmax,boundf,boundk,outfile,outputFreq);

  return ds->cvs[i];
}

restraint * DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr,
			 double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq )
{
  int i;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->rstr_count++;
  if (ds->rstr_count>ds->rstr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->rstr_count-1;
  ds->rstrs[i]=New_restraint(k,z,nCV,cvc,rftypstr,zmin,zmax,boundf,boundk,outfile, outputFreq);

  //restart
  if (outputFreq>0) {
    if (outputFreq < ds->rstr_save_frequency) {ds->rstr_save_frequency = outputFreq;}
  }

  return ds->rstrs[i];
}
