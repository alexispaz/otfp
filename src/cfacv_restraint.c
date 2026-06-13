#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfacv_internal.h"
#include "cfacv.h"
#include "measurements.h"

char * RFSTRINGS[NULL_RF] = {"HARMONIC","HARMCUTO"};
int rf_getityp ( char * typ ) {
  int i;
  for (i=0;i<NULL_RF&&strcmp(typ,RFSTRINGS[i]);i++);
  if (i<NULL_RF) return i;
  else return -1;
}

char * BADCVSTRINGMESSAGE = "NOT FOUND";
char * rf_getstyp ( int ityp ) {
  if (ityp<NULL_RF) return RFSTRINGS[ityp];
  else return BADCVSTRINGMESSAGE;
}

restraint * New_restraint (
    double k, double z, int nCV,
    double * cvc, char * rftypstr,
    double zmin, double zmax,
    char * boundstr, double boundk,
    char * outfile, int outputFreq
    ) {

  restraint * newr=calloc(1,sizeof(restraint));
  int i,b_pbc,b_peri;
  b_pbc=0;
  b_peri=0;

  newr->rfityp=rf_getityp(rftypstr);
  newr->k=k;
  newr->z=z;
  newr->u=newr->f=0.0;
  newr->nCV=nCV;
  newr->cvc=malloc(nCV*sizeof(double));
  if (cvc) for (i=0;i<nCV;i++) newr->cvc[i]=cvc[i];
  else for (i=0;i<nCV;i++) newr->cvc[i]=0;

  // evolve function
  newr->tamdOpt=NULL;
  newr->smdOpt=NULL;
  newr->adamOpt=NULL;

  newr->evolve=1;
  newr->evolveFunc = NULL;

  newr->tamd_noise=0.0;
  newr->tamd_restraint=0.0;

  // boundary function
  newr->min=zmin;
  newr->max=zmax;
  newr->half_domain=0.5*(zmax-zmin);
  newr->boundk=boundk;
  if     (!strcmp(boundstr,"SOFTUPPER")) {newr->boundFunc = SoftUpperWall      ;}
  else if(!strcmp(boundstr,"SOFTLOWER")) {newr->boundFunc = SoftLowerWall      ;}
  else if(!strcmp(boundstr,"SOFT"     )) {newr->boundFunc = SoftWalls          ;}
  else if(!strcmp(boundstr,"PBC"))       {newr->boundFunc = pbc; b_pbc=1       ;}
  else if(!strcmp(boundstr,"PERIODIC" )) {newr->boundFunc = Periodic; b_peri=1 ;}
  else if(!strcmp(boundstr,"NADA"     )) {newr->boundFunc = nada               ;}
  else {
    fprintf(stderr, "Error: boundary type not recognized");
    exit(1);
  }

  // energy function
  if     (b_pbc) {
    if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart_pbc;}
    else if (newr->rfityp==HARMCUTO) {newr->energyFunc = HarmonicCart_cutoff_pbc;}
  } else if (b_peri) {
    if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart_pbc;}
    else {
      fprintf(stderr, "Error: TODO HARMCUTO with periodic");
      exit(1);
    }
  } else {
    if      (newr->rfityp==HARMONIC) {newr->energyFunc = HarmonicCart;}
    else if (newr->rfityp==HARMCUTO) {newr->energyFunc = HarmonicCart_cutoff;}
  }

  // output
  newr->boutput=(outputFreq>0);
  if (newr->boutput) {
    newr->boutput=1;
    newr->outputFreq=outputFreq;
    newr->ofp=fopen(outfile,"w");
  }

  return newr;
}

void restraint_free ( restraint * r ) {
  if (!r) return;
  free(r->cvc);
  free(r->tamdOpt);
  free(r->smdOpt);
  free(r->adamOpt);
  if (r->boutput && r->ofp) fclose(r->ofp);
  free(r);
}

tamdOptStruct * New_tamdOptStruct ( double g, double kT, double dt, int riftyp) {
  tamdOptStruct * tamd=malloc(sizeof(tamdOptStruct));
  tamd->kT=kT;
  tamd->gamma=g;
  tamd->dt=dt;
  if (g>0.0) tamd->ginv=1.0/g;
  else tamd->ginv=0.0;
  tamd->noise = sqrt(2.0*kT*dt*tamd->ginv);
  return tamd;
}

smdOptStruct * New_smdOptStruct ( double target, int t0, int t1 ) {
  smdOptStruct * smd=malloc(sizeof(smdOptStruct));
  smd->t0=t0;
  smd->t1=t1;
  if (!(t1-t0)) smd->invinterval=0.;
  else smd->invinterval = 1.0/(t1-t0);
  smd->target=target;
  return smd;
}

adamOptStruct * New_adamOptStruct ( double a, double b1, double b2, double e, double decay, double kkT, double mvar, double vvar, double tmin) {
  adamOptStruct * adam=malloc(sizeof(adamOptStruct));

  adam->a=a;
  adam->decay=decay;

  adam->mvar=kkT*mvar;
  adam->vvar=kkT*kkT*vvar;
  adam->vmed=kkT;
  adam->zvar=a*a*mvar;

  adam->b1=b1;
  adam->b2=b2;
  adam->e=e;

  adam->t=1;
  adam->tmin=tmin;
  adam->stage=0;

  adam->v=0.;
  adam->m=0.;
  return adam;
}

int smdOptInit ( smdOptStruct * smd, double initval) {
  double dist=0.0;
  if (!smd) return -1;

  smd->initval=initval;
  dist = smd->target-smd->initval;
  smd->increment=dist*smd->invinterval;
  printf("OTFP: smd increment is %.5f based on target %.5f and initval %.5f\n",
         smd->increment,smd->target,smd->initval);
  return 0;
}

void ds_saverestrains ( DataSpace * ds, char * filename ) {
  int i;
  restraint * r;
  FILE * ofs;

  ofs=fopen(filename,"w");
  fwrite(&ds->K,sizeof(int),1,ofs);
  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    fwrite(&(r->z),sizeof(double),1,ofs);
    fwrite(&(r->val),sizeof(double),1,ofs);
  }

  fclose(ofs);
}

void ds_loadrestrains ( DataSpace * ds, char * filename ) {
  int i;
  restraint * r;
  FILE * ofs;

  ofs=fopen(filename,"r");

  fread(&i,sizeof(int),1,ofs);
  fprintf(stdout,"CFACV) Recovering restraint state of step %i of some previous simulation\n",i);
  if ( i!=ds->K ) {
    fprintf(stderr,"CFACV) ERROR restraint restart because holding object has not the proper size\n");
    exit(-1);
  }

  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    fread(&(r->z),sizeof(double),1,ofs);
    fread(&(r->val),sizeof(double),1,ofs);
    printf("OTFP: Recovering restraint %i, z=%.5f; val=%.5f \n",i,r->z,r->val);
  }

  fclose(ofs);
}

void restraint_reset_value ( restraint * r ) {
  r->val=0.0;
}

void restraint_add_cv_value ( restraint * r, int cv_index, double value ) {
  r->val+=r->cvc[cv_index]*value;
}

void restraint_init_target_from_value ( restraint * r ) {
  r->z=r->val;
  if (r->smdOpt) smdOptInit(r->smdOpt,r->val);
}

int restraint_compute_energy ( restraint * r ) {
  return r->energyFunc(r);
}

int restraint_should_output ( restraint * r, int timestep ) {
  return r->boutput && !(timestep % r->outputFreq);
}

int restraint_cv_count ( restraint * r ) {
  return r->nCV;
}

double restraint_cv_coefficient ( restraint * r, int cv_index ) {
  return r->cvc[cv_index];
}

double restraint_force ( restraint * r ) {
  return r->f;
}

double restraint_target ( restraint * r ) {
  return r->z;
}

int restraint_chapeau_id ( restraint * r ) {
  return r->chid;
}

int restraint_chapeau_dim ( restraint * r ) {
  return r->chdm;
}

int restraint_evolve_for_timestep ( restraint * r, int timestep ) {
  int done;

  r->boundFunc(r);

  if (r->tamdOpt) r->evolveFunc(r);

  if (r->smdOpt) {
    r->evolve=(int)(r->smdOpt->t0<=timestep)&&(r->smdOpt->t1>=timestep);
    r->evolveFunc(r);
  }

  if (r->adamOpt) {
    done=r->evolveFunc(r);
    if (done==1) return 1;
  }

  return 0;
}

void restraint_stop_evolution ( restraint * r ) {
  r->evolve=0;
}

int HarmonicCart ( restraint * r ) {
  double tmp=r->val-r->z;
  r->f=-r->k*tmp;
  r->u=0.5*r->k*tmp*tmp;
  return 0;
}

int HarmonicCart_cutoff ( restraint * r ) {
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  }

  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*gasdev();
    r->evolve=1;
  }

  HarmonicCart(r);
  return 0;
}

int HarmonicCart_pbc ( restraint * r ) {
  double tmp=r->val-r->z;
  if (tmp>r->half_domain) {tmp-=2*r->half_domain;}
  if (tmp<-r->half_domain) {tmp+=2*r->half_domain;}
  r->f=-r->k*tmp;
  r->u=0.5*r->k*tmp*tmp;
  return 0;
}

int HarmonicCart_cutoff_pbc ( restraint * r ) {
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  }

  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*gasdev();
    r->evolve=1;
  }

  HarmonicCart_pbc(r);
  return 0;
}

int nada ( restraint * r ) {
  return 0;
}

int SoftWalls ( restraint * r ) {
  SoftUpperWall(r);
  SoftLowerWall(r);
  return 0;
}

int SoftLowerWall ( restraint * r ) {
  double aux;
  aux=r->z-r->min;
  if (aux>0.) return 0;
  // r->f is minus the force on z!
  r->f+=r->boundk*aux;
  r->u+=.5*r->boundk*aux*aux;
  return 0;
}

int SoftUpperWall ( restraint * r ) {
  double aux;
  aux=r->z-r->max;
  if (aux<0.) return 0;
  // r->f is minus the force on z!
  r->f+=r->boundk*aux;
  r->u+=.5*r->boundk*aux*aux;
  return 0;
}

int Periodic ( restraint * r ) {
  if (r->z>M_PI) r->z-=2*M_PI;
  else if (r->z<-M_PI) r->z+=2*M_PI;
  return 0;
}

int pbc ( restraint * r ) {
  if (r->z>r->max) r->z-=2*r->half_domain;
  if (r->z<-r->min) r->z+=2*r->half_domain;
  return 0;
}

int cbd ( restraint * r ) {
  if (!r->evolve) return 0;

  tamdOptStruct * tamd=r->tamdOpt;
  double f = -r->f;
  double dd = tamd->ginv*tamd->dt*f;
  double rd = tamd->noise*gasdev();

  r->tamd_noise=rd;
  r->tamd_restraint=dd;
  r->z=r->z+dd+rd;
  return 0;
}

int adam ( restraint * r) {
  if (!r->evolve) return 1;

  adamOptStruct * adam=r->adamOpt;
  double e = adam->e;
  double b1 = adam->b1;
  double b2 = adam->b2;
  double dz;

  adam->t = adam->t + 1;
  adam->m=b1*adam->m+(1.-b1)*r->f;
  adam->v=b2*adam->v+(1.-b2)*r->f*r->f;

  if(adam->stage == 0) {
      b1 = pow(b1, adam->t);
      b2 = pow(b2, adam->t);
      double m;
      double v;

      m = adam->m / (1. - b1);
      v = adam->v / (1. - b2);

      dz=adam->a*m/(sqrt(v)+e);

      if(b2 < e) {
        fprintf(stdout,"ADAM) moving to stage 1 at step %i: ensure tmin=%i steps\n", adam->t,adam->tmin);
        adam->stage = 1;
      }

  } else if (adam->stage == 1){

      dz=adam->a*adam->m/(sqrt(adam->v)+e);

      if (adam->t>adam->tmin) {
        fprintf(stdout,"ADAM) moving to stage 2 at step %i: start termination check\n",adam->t);
        adam->stage = 2;
      }

  } else if (adam->stage == 2){

      dz=adam->a*adam->m/(sqrt(adam->v)+e);

      if (adam->m*adam->m < adam->mvar) {
         double aux=adam->v-adam->vmed;
         if (aux*aux < adam->vvar) {
            if (dz*dz < adam->zvar) {
              fprintf(stdout,"ADAM) moving to stage 3 at step %i: start decay learning rate\n",adam->t);
              adam->stage = 3;
              adam->t = 1;
            }
         }
      }

  } else {
    adam->a = adam->a/(1. + adam->decay * adam->t);
    dz=adam->a*adam->m/(sqrt(adam->v)+e);

    if(adam->a<e) {
      fprintf(stdout,"ADAM) terminating adam restrain after %i decay steps\n",adam->t);
      return 1;
    }
  }

  r->z=r->z-dz;

  return 0;
}

int uniformvelocity ( restraint * r ) {
  if (!r->evolve) return 0;

  double f=r->smdOpt->increment;
  r->z=r->z+f;
  return 0;
}

int restr_set_rchid ( restraint * r, DataSpace * ds, int chid) {
  if (!ds) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (!r) {fprintf(stderr,"CFACV) null argument\n"); exit(-1);}
  if (chid>=ds->ch_num) return -1;
  r->chid = chid;
  return 0;
}

int restr_AddTamdOpt ( restraint * r, double g, double kt, double dt, int chid, int chdm  ) {

  if (!r) return -1;

  r->evolveFunc = cbd;
  r->tamdOpt=New_tamdOptStruct(g,kt,dt,r->rfityp);
  r->chid = chid;
  r->chdm = chdm;
  return 0;
}

int restr_UpdateTamdOpt ( restraint * r, double g, double kt, double dt ) {
  tamdOptStruct * tamd;
  if (!r) return -1;
  tamd = r->tamdOpt;
  tamd->kT = kt;
  tamd->gamma = g;
  tamd->dt = dt;
  tamd->ginv = 1.0/g;
  tamd->noise = sqrt(2.0*kt*dt*tamd->ginv);
  return 0;
}

int restr_AddAdamOpt  ( restraint * r, double a, double b1, double b2, double e, double decay, double kT, double tmin) {
  if (!r) return -1;
  r->evolveFunc = adam;
  double aux=r->k*kT;
  r->adamOpt = New_adamOptStruct(a,b1,b2,e,decay,aux,(1.-b1)/(1.+b1),(1.-b2)/(1.+b2), tmin);
  return 0;
}

int restr_AddSmdOpt  ( restraint * r, double target, int t0, int t1 ) {
  if (!r) return -1;
  r->evolveFunc = uniformvelocity;
  r->smdOpt = New_smdOptStruct(target,t0,t1);
  r->chid = 0;
  return 0;
}

double restr_getz ( restraint * r ) {
  if (!r) return -1;
  return r->z;
}

double restr_getu ( restraint * r ) {
  if (!r) return -1;
  return r->u;
}

void restr_output ( restraint * r ) {
  fprintf(r->ofp,"%11.5f",r->z);
  fprintf(r->ofp," %11.5f",r->f);
  fprintf(r->ofp," %d\n",r->chid);
  fflush(r->ofp);
}
