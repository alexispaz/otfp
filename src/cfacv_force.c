#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "cfacv_internal.h"
#include "chapeau_internal.h"
#include "cvs_internal.h"

static void handler(int sig) {
  void *array[10];
  size_t size;

  size = backtrace(array, 10);

  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

static void DataSpace_clear_forces ( DataSpace * ds ) {
  int i,d;

  for (i=0;i<ds->ctr_total;i++) {
    for (d=0;d<3;d++) ds->ctr_data[i][d]=0.0;
  }
}

static void DataSpace_compute_restraint_values ( DataSpace * ds, int first ) {
  int i,j;
  int nCV;
  cv * cv;
  restraint * r;

  for (i=0;i<ds->rstr_total;i++) {
    r=ds->rstrs[i];
    nCV=restraint_cv_count(r);
    restraint_reset_value(r);
    for (j=0;j<nCV;j++) {
      cv=ds->cvs[j];
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

  for (i=0;i<ds->rstr_total;i++) {
    r=ds->rstrs[i];
    nCV=restraint_cv_count(r);

    restraint_compute_energy(r);
    if (restraint_should_output(r,timestep)) restr_output(r);

    for (j=0;j<nCV;j++) {
      cv=ds->cvs[j];
      coefficient=restraint_cv_coefficient(r,j);
      if (!coefficient) continue;

      for (jj=0;jj<cv_center_count(cv);jj++) {
        for (d=0;d<3;d++) {
          ds->ctr_data[cv_center_index(cv,jj)][d]+=coefficient*cv_gradient_component(cv,jj,d)*restraint_force(r);
        }
      }
    }
  }
}

static void DataSpace_apply_cv_boundary_forces ( DataSpace * ds, int timestep ) {
  int j,jj,d;
  cv * cv;

  for (j=0;j<ds->cv_total;j++) {
    cv=ds->cvs[j];

    cv_reset_force_energy(cv);
    cv_apply_boundary(cv);

    for (jj=0;jj<cv_center_count(cv);jj++) {
      for (d=0;d<3;d++) {
        ds->ctr_data[cv_center_index(cv,jj)][d]+=cv_gradient_component(cv,jj,d)*cv_force(cv);
      }
    }

    if (cv_should_output(cv,timestep)) cv_output(cv);
  }
}

static void DataSpace_update_fes ( DataSpace * ds, int timestep, double bias ) {
  int i;
  restraint * r;
  chapeau * ch;

  if (!ds->fes_enabled) return;

  if (timestep>ds->fes_begin_accum) {
    for (i=0;i<ds->rstr_total;i++) {
      r=ds->rstrs[i];
      ch=ds->fes[restraint_chapeau_id(r)];
      chapeau_set_sample(ch,restraint_chapeau_dim(r),restraint_target(r),-restraint_force(r)); // F=k(\theta-z)
    }
    for (i=0;i<ds->fes_total;i++) {
      chapeau_accumulate_sample(ds->fes[i],bias);
    }
  }

  if (!ds->fes_evolve_parameters) return;

  for (i=0;i<ds->fes_total;i++) {
    ch=ds->fes[i];

    if (chapeau_should_solve(ch,timestep,ds->fes_begin_solve)) {
      chapeau_solve(ch);
    }

    if (chapeau_should_output(ch,timestep)) {
      chapeau_output(ch,timestep);

      // TODO: separete savestate from outputFreq
      chapeau_savestate_default(ch);
    }
  }
}

static int DataSpace_evolve_rstrs ( DataSpace * ds, int timestep ) {
  int i;
  int rabort=0;
  restraint * r;

  for (i=0;i<ds->rstr_total;i++) {
    r=ds->rstrs[i];

    if (restraint_evolve_for_timestep(r,timestep)) {
      ds_saverestrains(ds,ds->filename);
      restraint_stop_evolution(r);
      rabort++;
    }
  }

  if(ds->rstr_total>0) {
    if (!(timestep % ds->rstr_save_frequency)) ds_saverestrains(ds,ds->filename);
  }

  if(rabort==ds->rstr_total) return 2;
  return 0;
}

int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep, double bias ) {
  signal(SIGINT, handler);
  if (!ds) return -1;

  DataSpace_clear_forces(ds);
  DataSpace_compute_restraint_values(ds,first);
  DataSpace_apply_restraint_forces(ds,timestep);
  DataSpace_apply_cv_boundary_forces(ds,timestep);
  DataSpace_update_fes(ds,timestep,bias);

  return DataSpace_evolve_rstrs(ds,timestep);
}
