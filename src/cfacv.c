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
    cv_compute(ds->cv[i],ds->R);
    bias+=cv_bias_value(ds->cv[i]);
  }
 
  return bias;
}

static void DataSpace_clear_forces ( DataSpace * ds ) {
  int i,d;

  for (i=0;i<ds->N;i++) {
    for (d=0;d<3;d++) ds->R[i][d]=0.0;
  }
}

static void DataSpace_compute_restraint_values ( DataSpace * ds, int first ) {
  int i,j;
  int nCV;
  cv * cv;
  restraint * r;

  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    nCV=restraint_cv_count(r);
    restraint_reset_value(r);
    for (j=0;j<nCV;j++) {
      cv=ds->cv[j];
      restraint_add_cv_value(r,j,cv_value(cv));
    }

    if (first) restraint_init_target_from_value(r);
  }
}

static void DataSpace_apply_restraint_forces ( DataSpace * ds, int timestep ) {
  int i,j,jj,d;
  int nCV;
  double coefficient;
  cv * cv;
  restraint * r;

  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];
    nCV=restraint_cv_count(r);

    restraint_compute_energy(r);
    if (restraint_should_output(r,timestep)) restr_output(r);

    for (j=0;j<nCV;j++) {
      cv=ds->cv[j];
      coefficient=restraint_cv_coefficient(r,j);
      if (!coefficient) continue;

      for (jj=0;jj<cv_center_count(cv);jj++) {
        for (d=0;d<3;d++) {
          ds->R[cv_center_index(cv,jj)][d]+=coefficient*cv_gradient_component(cv,jj,d)*restraint_force(r);
        }
      }
    }
  }
}

static void DataSpace_apply_cv_boundary_forces ( DataSpace * ds, int timestep ) {
  int j,jj,d;
  cv * cv;

  for (j=0;j<ds->ncv;j++) {
    cv=ds->cv[j];

    cv_reset_force_energy(cv);
    cv_apply_boundary(cv);

    for (jj=0;jj<cv_center_count(cv);jj++) {
      for (d=0;d<3;d++) {
        ds->R[cv_center_index(cv,jj)][d]+=cv_gradient_component(cv,jj,d)*cv_force(cv);
      }
    }

    if (cv_should_output(cv,timestep)) cv_output(cv);
  }
}

static void DataSpace_update_chapeau ( DataSpace * ds, int timestep, double bias ) {
  int i;
  restraint * r;
  chapeau * ch;

  if (!ds->doAnalyticalCalc) return;

  if (timestep>ds->beginaccum) {
    for (i=0;i<ds->K;i++) {
      r=ds->restr[i];
      ch=ds->ch[restraint_chapeau_id(r)];
      chapeau_set_sample(ch,restraint_chapeau_dim(r),restraint_target(r),-restraint_force(r)); // F=k(\theta-z)
    }
    for (i=0;i<ds->ch_num;i++) {
      chapeau_accumulate_sample(ds->ch[i],bias);
    }
  }

  if (!ds->evolveAnalyticalParameters) return;

  for (i=0;i<ds->ch_num;i++) {
    ch=ds->ch[i];

    if (chapeau_should_solve(ch,timestep,ds->beginsolve)) {
      chapeau_solve(ch);
    }

    if (chapeau_should_output(ch,timestep)) {
      chapeau_output(ch,timestep);

      // TODO: separete savestate from outputFreq
      chapeau_savestate_default(ch);
    }
  }
}

static int DataSpace_evolve_restraints ( DataSpace * ds, int timestep ) {
  int i;
  int rabort=0;
  restraint * r;

  for (i=0;i<ds->K;i++) {
    r=ds->restr[i];

    if (restraint_evolve_for_timestep(r,timestep)) {
      ds_saverestrains(ds,ds->filename);
      restraint_stop_evolution(r);
      rabort++;
    }
  }

  if(ds->K>0) {
    if (!(timestep % ds->restrsavefreq)) ds_saverestrains(ds,ds->filename);
  }

  if(rabort==ds->K) return 2;
  return 0;
}

int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep, double bias ) {
  // install SIGINT handler
  signal(SIGINT, handler);
  if (!ds) return -1;

  DataSpace_clear_forces(ds);
  DataSpace_compute_restraint_values(ds,first);
  DataSpace_apply_restraint_forces(ds,timestep);
  DataSpace_apply_cv_boundary_forces(ds,timestep);
  DataSpace_update_chapeau(ds,timestep,bias);

  return DataSpace_evolve_restraints(ds,timestep);
}

double DataSpace_RestraintEnergy ( DataSpace * ds ) {
  if (ds) {
    double u=0.0;
    int i;
    restraint * r;
    for (i=0;i<ds->K;i++) {
      r=ds->restr[i];
      u+=restr_getu(r);
    }
    return u;
  }
  return 0.0;    
}

void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res ) {
  if (ds) {
    int i;
    for (i=0;i<ds->ncv;i++) {
      if (active[i]) res[i]=cv_value(ds->cv[i]);
      else res[i]=0.0;
    }
  }
}

void DataSpace_ReportAll ( DataSpace * ds ) {
  if (ds) {
    int i;
    printf("CFACV_OTFP/C) DataSpace: %i Centers, %i CV's, %i Restraints\n",
	   ds->N, ds->ncv, ds->K);
    for (i=0;i<ds->N;i++) {
      printf("CFACV_OTFP/C)  Center %i : %.5f %.5f %.5f\n",
	     i,ds->R[i][0],ds->R[i][1],ds->R[i][2]);
    }
    for (i=0;i<ds->ncv;i++) {
      cv_report(ds->cv[i],i);
    }
    for (i=0;i<ds->K;i++) {
      restraint_report(ds->restr[i],i);
    }
    
  }
}
