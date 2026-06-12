/* Collective Variables module for NAMD via tclForces
 * Cameron F. Abrams
 *
 * Extended for 2D OTFP and RE OTFP
 * Sergio A. Paz
 *
 * This header declares structured data types and functions for
 * handling collective variables and their restraints. 
 * 
 * The central concept is that of a "restraint" which is applied to
 * one or a linear combination of "collective variables".  A restraint
 * computes a force on a collective variable ($\theta$) based on its
 * value and the value of a "target" or "auxiliary" variable ($z$).
 * The auxiliary variable can be constant, or driven according to
 * steered MD or temperature-accelerated MD.  It is the job of the
 * cfacv module to compute these forces after receiving from the NAMD
 * tclforces "calcforces" function the positions of all atoms or
 * atom-groups that participate in each CV owned by each restraint.
 * It is the job of the tclforces module to receive those forces and
 * transmit them to the atoms. This two-way communication happens
 * once per timestep on the master node during NAMD execution.
 *
 * In this particular OTFP impelementation, we are using linear chapeau
 * functions as the basis functions. We use this to compute effective pair
 * potentials between centers of mass of atomgroups (like water molecules).
 *
 */
 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cfacv_obj.h"
#include "chapeau.h"       // basis functions
#include "measurements.h"  // random numbers
 
#include "cvs.h"
#include "cfacv.h"
                    
// Handling SIGSEGV or SIGINT signal to make a call stack using backtrace.
// The call stack can be used to profile the program
// See for instance:
// http://stackoverflow.com/a/378024/1342186
// http://stackoverflow.com/a/77336/1342186

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void cfacvBanner ( void ) {
  printf("OTFP: Version (Git Hash): %s\n",VERSION);
}

FILE * my_fopen ( char * name, char * code ) {
  if (strcmp(name,"stdout")) return fopen(name,code);
  else return stdout;
}

FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds ) {
  FILE * fp=fopen(name,code);
  fwrite(&outputLevel,sizeof(unsigned int),1,fp);
  fwrite(&(ds->K),sizeof(int),1,fp);
  return fp;
}


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




// RESTRAIN STRUCTURES

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
  //FIXME: Inquire if the restrain is linked to periodic instead.
  //if (periodic) {
  //  if (dist<-M_PI) {
  //    dist+=2*M_PI;
  //  }
  //  else if (dist > M_PI) {
  //    dist-=2*M_PI;
  //  } 
  //}
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

//ENERGY FUNCTIONS

int HarmonicCart ( restraint * r ) {
  double tmp=r->val-r->z;
  r->f=-r->k*tmp;
  r->u=0.5*r->k*tmp*tmp;
  return 0;
}

int HarmonicCart_cutoff ( restraint * r ) {
  //double r1=10.0,r2=12.0;
 
  // Non tie
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  } 


  //// Set up initial position above the CV
  //if(!r->evolve) {
  //  r->z=r->val;
  //  r->evolve=1;
  //  r->f=0.0;
  //  r->u=0.0;
  //  return 0;
  //}

  // Set up initial acording to a distribution around the CV
  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*gasdev();
    r->evolve=1;
  }

  // Abrupt cutoff
  HarmonicCart(r);

  //// Force an energy with a smooth cutoff
  //if(r->z<r1) {
  //  aux=0.0;
  //  aux1=1.0;
  //} else {
  //tmp==r->val-r->z;
  //  if (r->pbc) {
  //    if (tmp>r->half_domain) tmp-=2*r->half_domain;
  //    if (tmp<-r->half_domain) tmp+=2*r->half_domain;
  //  }
  //  aux=M_PI*(tmp-(r1+r2)*0.5)/(r2-r1);
  //  aux1=0.5-0.5*sin(aux);
  //  aux=-cos(aux)*M_PI*tmp/(r2-r1)*0.5;
  //  //dfcut=-0.5_dp*cos(aux)*pio2*tmp/(r2-r1);
  //}
  //aux2=0.5*r->k*tmp*tmp;
  //r->f=-r->k*tmp*aux1-aux2*aux;
  //r->u=aux2*aux1;

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

  // Release tie
  if (r->z>r->max) {
    r->f=0.0;
    r->u=0.0;
    r->evolve=0;
    return 0;
  } 

  // Set up initial acording to a distribution around the CV
  if(!r->evolve) {
    r->z=r->val+sqrt(r->tamdOpt->kT/r->k)*gasdev();
    r->evolve=1;
  }

  // Abrupt cutoff
  HarmonicCart_pbc(r);

  return 0;
}
      



//BOUNDARIES
 
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
  // TODO: Periodic is the same that PBC and can be unified.
  if (r->z>M_PI) r->z-=2*M_PI;
  else if (r->z<-M_PI) r->z+=2*M_PI;
  return 0;
}
 
int pbc ( restraint * r ) {
  if (r->z>r->max) r->z-=2*r->half_domain;
  if (r->z<-r->min) r->z+=2*r->half_domain;
  return 0;
}


// EVOLUTION OF RESTRIANS

/* Brownian dynamics (TAMD)*/
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
       
/* Adam (adam)*/
// Kingma, D.P., Ba, J., 2017. Adam: A Method for
// Stochastic Optimization.
// https://doi.org/10.48550/arXiv.1412.6980
int adam ( restraint * r) {
  if (!r->evolve) return 1;

  adamOptStruct * adam=r->adamOpt;
  double e = adam->e;
  double b1 = adam->b1;
  double b2 = adam->b2;
  double dz;

  // timestep
  adam->t = adam->t + 1;
      
  //Update biased first moment estimate
  adam->m=b1*adam->m+(1.-b1)*r->f;

  //Update biased second raw moment estimate
  adam->v=b2*adam->v+(1.-b2)*r->f*r->f;

  // Compute bias-corrections. Avoid if possible, considering that for long t,
  // correction tends to one.

  if(adam->stage == 0) {
      // Calculate the maximum of b1^t and b2^t
      b1 = pow(b1, adam->t);
      b2 = pow(b2, adam->t);
      double m;
      double v;

      // Bias corrections
      m = adam->m / (1. - b1);
      v = adam->v / (1. - b2);
      
      dz=adam->a*m/(sqrt(v)+e);

      // b2 is always bigger than b1
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

      // Teste for Termination.
      // Termination implies to start the learning rate decay.
      // In the minima, is probable that fluctuations are due to noise, and not a slope.
      // Assuming the minima behave as a SHO of constant k: 
      //    Var(g) = E(g^2)-E(g)^2 = k*k_BT - 0 = k*k_BT
      // On the other hand, as in eq 1 of Kingma2017:
      //    m=(1-b_1)\sum_i=1^t b_1^(t-i)g_i
      // Since g_i are independent, variance are the sum of variance:
      //   Var(m)=(1-b_1)^2\sum_i=1^t b_1^{2(t-i)}Var(g_i)
      //         =(1-b_1)^2 k*k_BT \sum_i=1^t b_1^{2(t-i)}
      // Using the geometrical serie
      //         =(1-b_1)^2 k*k_BT (1-b_1^{2t})/(1-b_1^2)
      // Taking b_1^{2t}->0 for long t
      //         = k*k_BT (1-b_1)/(1+b_1)
      // I use those equations to build one of the termination criteria. 
      if (adam->m*adam->m < adam->mvar) {
         
         // Do the same with v
         double aux=adam->v-adam->vmed;
         if (aux*aux < adam->vvar) {
         
            // Do the same with z
            if (dz*dz < adam->zvar) {
              fprintf(stdout,"ADAM) moving to stage 3 at step %i: start decay learning rate\n",adam->t);
              adam->stage = 3;
              adam->t = 1;
            }
         }
      }

  } else {
    // Adding a learning rate decay
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
     
/* Constant Velocity (SMD)*/
int uniformvelocity ( restraint * r ) {
  if (!r->evolve) return 0;

  double f=r->smdOpt->increment;
  r->z=r->z+f;
  return 0;
}

//double smdUpdate ( double z, double increment, int OK ) {
//  return z+OK*increment;
//}
//
//double smdUpdate_Periodic ( double z, double increment, int OK ) {
//  double newz=z+increment;
//  if (newz>M_PI) newz-=2*OK*M_PI;
//  else if (newz<-M_PI) newz+=2*OK*M_PI;
//  return newz;
//}







// DATA SPACE


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
  // rescalevels [expr sqrt(1.0*$NEWTEMP/$OLDTEMP)] no needed in cbd
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
  // fprintf(r->ofp,"%11.5f",r->val);
  fprintf(r->ofp," %11.5f",r->f);
  fprintf(r->ofp," %d\n",r->chid);
  fflush(r->ofp);
}

double DataSpace_ComputeCVs ( DataSpace * ds ) {
  int i;
  double bias;
 
  bias=0;
  if (!ds) return -1;
  for (i=0;i<ds->ncv;i++){
    ds->cv[i]->calc(ds->cv[i],ds->R);
    if (ds->cv[i]->amd) bias+=ds->cv[i]->val;
  }
 
  return bias;
}


int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep, double bias ) {
  int i,j,jj;
  double * cvc;
  int d=0;
  int rabort;
  int K=ds->K;
  cv * cv;
  restraint * r;
  chapeau * ch;
      
  // install SIGINT handler
  signal(SIGINT, handler);   
  if (!ds) return -1;

  /* clear out the position arrays to hold forces to be communicated
     back to tclForces */
  for (i=0;i<ds->N;i++) for (d=0;d<3;d++) ds->R[i][d]=0.0;

  /* Loop over all restraints and compute restraint values from their CV's */
  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    cvc=r->cvc;
    r->val=0.0;
    for (j=0;j<r->nCV;j++) {
      cv=ds->cv[j];
      r->val+=cvc[j]*cv->val;
    }
    
    if (first) {
      r->z=r->val;
      if (r->smdOpt) smdOptInit(r->smdOpt,r->val);
    }
  }

  //printf("OTFP: %i restraint %i val %.4f targ %.4f\n",timestep,i,r->val,r->z);

  /* compute forces on all restraints */
  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    cvc=r->cvc;

    // Compute the energy and the force of the restrain
    r->energyFunc(r);

    // restraint output
    if (r->boutput) {
      if (!(timestep % r->outputFreq)) restr_output(r);
    }

    /* accumulate force increments in ds->R for each real particle*/
    /* cv->gr[jj][kk] = \partial CV / \partial (kk-coord of center jj of CV) */

    // For each CV in this restraint
    for (j=0;j<r->nCV;j++) {
      cv=ds->cv[j];
      if (r->cvc[j]) {
        
        // For each center in this CV
        for (jj=0;jj<cv->nC;jj++) {

          /* Increment the forces of the jj'th cartesian center on the j'th cv of this restraint */
          for (d=0;d<3;d++) ds->R[ cv->ind[jj] ][d]+=r->cvc[j]*cv->gr[jj][d]*r->f;
        }
      }
    }
  }

  // For each CV in dataspace
  for (j=0;j<ds->ncv;j++) {
    cv=ds->cv[j];

    // Restart forces in CVs
    cv->f=0;
    cv->u=0;

    // Compute the boundary forces
    cv->boundFunc(cv);
    
    // Add boundary forces
    for (jj=0;jj<cv->nC;jj++) {
      for (d=0;d<3;d++) ds->R[ cv->ind[jj] ][d]+=cv->gr[jj][d]*cv->f;
    }

    // cv output
    if (cv->boutput) {
      if (!(timestep % cv->outputFreq)) cv_output(cv);
    }
     
  }
   
  // OTFP
  if (ds->doAnalyticalCalc) {
 
    // Add statistic to b and A matrix
    if (timestep>ds->beginaccum) {
      for (i=0;i<ds->K;i++) {
        r=ds->restr[i];
        ch=ds->ch[r->chid];
        ch->r[r->chdm]=r->z;
        ch->f[r->chdm]=-r->f; // F=k(\theta-z)
      }
      for (i=0;i<ds->ch_num;i++) {
        ds->ch[i]->accumulate(ds->ch[i],bias);
      }
    }

    // Solving FES
    if (ds->evolveAnalyticalParameters) {

      for (i=0;i<ds->ch_num;i++) {
        ch=ds->ch[i];

        if (timestep>ds->beginsolve && !(timestep % ch->nupdate)) {
          // chapeau_savestate(ch,"after.ch0\0");
         // getchar();
          chapeau_solve(ch);

          //for (j=0;j<K;j++) {
          //  r=ds->restr[j];
          //  r->ofp=freopen(NULL, "w", r->ofp);
          //}
          // chapeau_savestate(ch,"before.ch0\0");
        }

        //aux=ch->nsample/ch->nupdate
        //if(!(aux%ch->outputFreq)) chapeau_output(ch);

        if (!(timestep % ch->outputFreq)) {
          chapeau_output(ch,timestep);

          // TODO: separete savestate from outputFreq
          chapeau_savestate(ch,ch->filename);
        }
      }
    }
  }
   
  // Evolving auxilary variables
  rabort=0;
  for (i=0;i<K;i++) {
    r=ds->restr[i];

    // Compute the boundary forces not included in the atom forces!!
    r->boundFunc(r);

    // Evolution of TAMD restraints
    if (r->tamdOpt) r->evolveFunc(r);
  
    // Evolution of SMD restraints
    if (r->smdOpt) {
      r->evolve=(int)(r->smdOpt->t0<=timestep)&&(r->smdOpt->t1>=timestep);
      r->evolveFunc(r);
    }
    if (r->adamOpt) {
      j=r->evolveFunc(r);
      if(j==1) {
        ds_saverestrains(ds,ds->filename);
        r->evolve=0;
        rabort=rabort+1;
      }
    }
  }

  //Save restr trajectory for future restart
  if(K>0) {
    if (!(timestep % ds->restrsavefreq)) ds_saverestrains(ds,ds->filename);
  }
 
  // Abort simulation when all restraint want it
  // this is usefull for things like adam.
  // To skip, add a dummy restraint without adam.
  if(rabort==K) {return 2;}
              
  return 0;
}

double DataSpace_RestraintEnergy ( DataSpace * ds ) {
  if (ds) {
    double u=0.0;
    int i;
    restraint * r;
    for (i=0;i<ds->K;i++) {
      r=ds->restr[i];
      u+=r->u;
    }
    return u;
  }
  return 0.0;    
}

void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res ) {
  if (ds) {
    int i;
    for (i=0;i<ds->ncv;i++) {
      if (active[i]) res[i]=ds->cv[i]->val;
      else res[i]=0.0;
    }
  }
}

void DataSpace_ReportAll ( DataSpace * ds ) {
  if (ds) {
    int i,j,k;
    printf("CFACV_OTFP/C) DataSpace: %i Centers, %i CV's, %i Restraints\n",
	   ds->N, ds->ncv, ds->K);
    for (i=0;i<ds->N;i++) {
      printf("CFACV_OTFP/C)  Center %i : %.5f %.5f %.5f\n",
	     i,ds->R[i][0],ds->R[i][1],ds->R[i][2]);
    }
    for (i=0;i<ds->ncv;i++) {
      printf("CFACV_OTFP/C)      CV %i : typ %s val %.5f ind ",i,cv_getstyp(ds->cv[i]->typ),
	     ds->cv[i]->val);
      for (j=0;j<ds->cv[i]->nC;j++) printf("%i ",ds->cv[i]->ind[j]);
      printf(" grads: ");
      for (j=0;j<ds->cv[i]->nC;j++) 
	for (k=0;k<3;k++) printf("%.5f ",ds->cv[i]->gr[j][k]);
      printf("\n");
    }
    for (i=0;i<ds->K;i++) {
      printf("CFACV_OTFP/C)       R %i : k %.3f z %.3f cvc ",i,ds->restr[i]->k,ds->restr[i]->z);
      for (j=0;j<ds->restr[i]->nCV;j++) printf("%.2f ",ds->restr[i]->cvc[j]);
      if (ds->restr[i]->tamdOpt) {
	printf("TAMD(kT=%.2f,gamma=%.2f,dt=%.3f) ",
	       ds->restr[i]->tamdOpt->kT,
	       ds->restr[i]->tamdOpt->gamma,
	       ds->restr[i]->tamdOpt->dt);
      }
      printf(" typ %s min %.5f max %.5f\n",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max);
      printf("\n");
    }
    
  }
}
