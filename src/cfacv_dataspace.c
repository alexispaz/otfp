#include <stdio.h>
#include <stdlib.h>

#include "cfacv_internal.h"
#include "cfacv.h"
#include "chapeau_internal.h"
#include "cvs.h"
#include "measurements.h"

static void OtfpStepState_init ( OtfpStepState * step, int nctr ) {
  int i;

  step->ctr_count=nctr;
  step->ctr_positions=malloc(nctr*sizeof(double*));
  step->ctr_forces=malloc(nctr*sizeof(double*));
  for (i=0;i<nctr;i++) {
    step->ctr_positions[i]=calloc(3,sizeof(double));
    step->ctr_forces[i]=calloc(3,sizeof(double));
  }
}

static void OtfpStepState_free ( OtfpStepState * step ) {
  int i;

  if (step->ctr_positions) {
    for (i=0;i<step->ctr_count;i++) free(step->ctr_positions[i]);
    free(step->ctr_positions);
  }
  if (step->ctr_forces) {
    for (i=0;i<step->ctr_count;i++) free(step->ctr_forces[i]);
    free(step->ctr_forces);
  }
}

static double * OtfpStepState_position ( OtfpStepState * step, int ctr ) {
  return step->ctr_positions[ctr];
}

static double * OtfpStepState_force ( OtfpStepState * step, int ctr ) {
  return step->ctr_forces[ctr];
}

static void OtfpStepState_clear_forces ( OtfpStepState * step ) {
  int i,d;

  for (i=0;i<step->ctr_count;i++) {
    for (d=0;d<3;d++) step->ctr_forces[i][d]=0.0;
  }
}

static void OtfpStepState_add_force ( OtfpStepState * step, int ctr, int dim, double force ) {
  step->ctr_forces[ctr][dim]+=force;
}

static void OtfpRegistry_init ( OtfpRegistry * registry, int nctr, int ncv, int nrstr ) {
  registry->ctr_total=nctr;
  registry->cv_total=ncv;
  registry->rstr_total=nrstr;
  registry->ctr_count=0;
  registry->cv_count=0;
  registry->rstr_count=0;

  registry->ctrs=calloc(nctr,sizeof(atomCenterStruct*));
  registry->cvs=calloc(ncv,sizeof(cv*));
  registry->rstrs=calloc(nrstr,sizeof(restraint*));
}

static void OtfpRegistry_free ( OtfpRegistry * registry ) {
  int i;

  if (registry->ctrs) {
    for (i=0;i<registry->ctr_total;i++) atomCenterStruct_free(registry->ctrs[i]);
    free(registry->ctrs);
  }
  if (registry->cvs) {
    for (i=0;i<registry->cv_total;i++) cv_free(registry->cvs[i]);
    free(registry->cvs);
  }
  if (registry->rstrs) {
    for (i=0;i<registry->rstr_total;i++) restraint_free(registry->rstrs[i]);
    free(registry->rstrs);
  }
}

static void OtfpFesState_init ( OtfpFesState * fes ) {
  fes->total=0;
  fes->enabled=0;
  fes->evolve_parameters=0;
  fes->begin_accum=0;
  fes->begin_solve=0;
  fes->chapeaux=NULL;
}

static void OtfpFesState_free ( OtfpFesState * fes ) {
  int i;

  if (fes->chapeaux) {
    for (i=0;i<fes->total;i++) chapeau_free(fes->chapeaux[i]);
    free(fes->chapeaux);
  }
}

DataSpace * NewDataSpace ( int Nctr, int Ncv, int Nrstr, long int seed ) {
  DataSpace * ds = calloc(1,sizeof(DataSpace));

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

  OtfpRegistry_init(&ds->registry,Nctr,Ncv,Nrstr);
  OtfpStepState_init(&ds->step,Nctr);
  OtfpFesState_init(&ds->fes);

  sprintf(ds->filename, "checkpoint.rstr");
  ds->rstr_save_frequency = 100000;

  return ds;
}

void DataSpace_free ( DataSpace * ds ) {
  if (!ds) return;
  OtfpStepState_free(&ds->step);
  OtfpRegistry_free(&ds->registry);
  OtfpFesState_free(&ds->fes);
  free(ds);
}

chapeau * DataSpace_get_chapeau ( DataSpace * ds, int i ) {
  return ds->fes.chapeaux[i];
}

int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double*  min, int * nKnots,
    double * max, int * periodic, int fes_begin_accum, int fes_begin_solve, char * outfile,
    int outfreq, int outlevel, int nupdate) {
  // Initialize FES state. The current implementation stores each FES state in
  // a chapeau object because chapeau is the active basis/kernel.
  int i;
  char buf[80], buf2[80];

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->fes.enabled=1;

  // This is to let the system equilibrate
  ds->fes.evolve_parameters=1;
  ds->fes.begin_accum=fes_begin_accum;
  ds->fes.begin_solve=fes_begin_solve;
  if (fes_begin_accum < 0) ds->fes.evolve_parameters=0;

  ds->fes.total=numrep;
  ds->fes.chapeaux=malloc(ds->fes.total*sizeof(chapeau*));

  for (i=0;i<ds->fes.total;i++) {
    ds->fes.chapeaux[i]=chapeau_alloc(dm,min,max,nKnots,periodic);
    sprintf(buf, "%s_ch%d.bsp", outfile,i);
    sprintf(buf2, "%s.ch%d", outfile,i);
    chapeau_setupoutput(ds->fes.chapeaux[i],buf,buf2,outfreq,outlevel);
    chapeau_set_update_frequency(ds->fes.chapeaux[i],nupdate);
  }

  return 0;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  return ds->registry.ctr_total;
}

int DataSpace_ctr_count ( DataSpace * ds ) {
  return ds->registry.ctr_total;
}

int DataSpace_cv_count ( DataSpace * ds ) {
  return ds->registry.cv_total;
}

int DataSpace_rstr_count ( DataSpace * ds ) {
  return ds->registry.rstr_total;
}

int DataSpace_fes_count ( DataSpace * ds ) {
  return ds->fes.total;
}

int DataSpace_fes_is_enabled ( DataSpace * ds ) {
  return ds->fes.enabled;
}

int DataSpace_fes_evolves ( DataSpace * ds ) {
  return ds->fes.evolve_parameters;
}

int DataSpace_fes_accum_begin ( DataSpace * ds ) {
  return ds->fes.begin_accum;
}

int DataSpace_fes_solve_begin ( DataSpace * ds ) {
  return ds->fes.begin_solve;
}

int DataSpace_rstr_save_frequency ( DataSpace * ds ) {
  return ds->rstr_save_frequency;
}

char * DataSpace_rstr_checkpoint_file ( DataSpace * ds ) {
  return ds->filename;
}

cv * DataSpace_cv ( DataSpace * ds, int index ) {
  return ds->registry.cvs[index];
}

restraint * DataSpace_rstr ( DataSpace * ds, int index ) {
  return ds->registry.rstrs[index];
}

chapeau * DataSpace_fes ( DataSpace * ds, int index ) {
  return ds->fes.chapeaux[index];
}

double ** DataSpace_ctr_positions ( DataSpace * ds ) {
  return ds->step.ctr_positions;
}

double * DataSpace_centerPosition ( DataSpace * ds, int i ) {
  return DataSpace_ctr_position(ds,i);
}

double * DataSpace_centerForce ( DataSpace * ds, int i ) {
  return DataSpace_ctr_force(ds,i);
}

double * DataSpace_ctr_position ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->registry.ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return OtfpStepState_position(&ds->step,i);
}

double * DataSpace_ctr_force ( DataSpace * ds, int i ) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (i>=ds->registry.ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}
  return OtfpStepState_force(&ds->step,i);
}

void DataSpace_clear_ctr_forces ( DataSpace * ds ) {
  OtfpStepState_clear_forces(&ds->step);
}

void DataSpace_add_ctr_force ( DataSpace * ds, int ctr, int dim, double force ) {
  OtfpStepState_add_force(&ds->step,ctr,dim,force);
}

int DataSpace_checkdata ( DataSpace * ds ) {
  int i,j;

  if (!ds) return 1;

  for (i=0;i<ds->registry.ctr_total;i++) {
    for (j=0;j<3;j++) {
      if (ds->step.ctr_positions[i][j]!=ds->step.ctr_positions[i][j]) {
        return 1;
      }
      if (ds->step.ctr_forces[i][j]!=ds->step.ctr_forces[i][j]) {
        return 1;
      }
    }
  }

  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->registry.ctr_total;i++) {
      fprintf(stdout,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->step.ctr_positions[i][j]);
      fprintf(stdout,"| ");
      for (j=0;j<3;j++)
	fprintf(stdout,"%.5f ",ds->step.ctr_forces[i][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
  }
  return 0;
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  int i,j;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->registry.ctr_count++;
  if (ds->registry.ctr_count>ds->registry.ctr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->registry.ctr_count-1;
  ds->registry.ctrs[i]=New_atomCenterStruct(n);
  for (j=0;j<n;j++) {
    ds->registry.ctrs[i]->ind[j]=ind[j];
    ds->registry.ctrs[i]->m[j]=m[j];
    ds->registry.ctrs[i]->M+=m[j];
  }

  return i;
}

cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
		      double zmin, double zmax,char * boundf,double boundk,
          char * outfile, int outputFreq )  {
  int i;
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->registry.cv_count++;
  if (ds->registry.cv_count>ds->registry.cv_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->registry.cv_count-1;
  ds->registry.cvs[i]=cv_init(typ,nind,ind,zmin,zmax,boundf,boundk,outfile,outputFreq);

  return ds->registry.cvs[i];
}

restraint * DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr,
			 double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq )
{
  int i;

  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}

  ds->registry.rstr_count++;
  if (ds->registry.rstr_count>ds->registry.rstr_total) {fprintf(stderr,"CFACV) out of size\n"); exit(-1);}

  i=ds->registry.rstr_count-1;
  ds->registry.rstrs[i]=New_restraint(k,z,nCV,cvc,rftypstr,zmin,zmax,boundf,boundk,outfile, outputFreq);

  //restart
  if (outputFreq>0) {
    if (outputFreq < ds->rstr_save_frequency) {ds->rstr_save_frequency = outputFreq;}
  }

  return ds->registry.rstrs[i];
}
