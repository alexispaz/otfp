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

#include "cfacv_internal.h"
#include "chapeau_internal.h"       // basis functions
#include "measurements.h"  // random numbers
 
#include "cvs_internal.h"
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
