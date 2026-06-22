#include <stdlib.h>

#include "cfacv.h"
#include "otfp_engine.h"

struct otfp_engine {
  DataSpace * ds;
  restraint ** restraints;
  int nrestr;
  int irestr;
};

otfp_engine * otfp_engine_create ( int ncenters, int ncvs, int nrestr, long seed ) {
  otfp_engine * engine = calloc(1,sizeof(otfp_engine));

  if (!engine) return NULL;

  engine->ds=NewDataSpace(ncenters,ncvs,nrestr,seed);
  if (!engine->ds) {
    free(engine);
    return NULL;
  }

  if (nrestr>0) {
    engine->restraints=calloc(nrestr,sizeof(restraint *));
    if (!engine->restraints) {
      DataSpace_free(engine->ds);
      free(engine);
      return NULL;
    }
  }
  engine->nrestr=nrestr;

  return engine;
}

void otfp_engine_free ( otfp_engine * engine ) {
  if (!engine) return;
  DataSpace_free(engine->ds);
  free(engine->restraints);
  free(engine);
}

int otfp_engine_add_center ( otfp_engine * engine, int natoms, const int * indices, const double * masses ) {
  if (!engine || !engine->ds || !indices || !masses) return -1;
  return DataSpace_AddAtomCenter(engine->ds,natoms,(int *)indices,(double *)masses);
}

int otfp_engine_set_center_position ( otfp_engine * engine, int center, const double xyz[3] ) {
  double * pos;

  if (!engine || !engine->ds || !xyz) return -1;

  pos=DataSpace_centerPosition(engine->ds,center);
  pos[0]=xyz[0];
  pos[1]=xyz[1];
  pos[2]=xyz[2];

  return 0;
}

const double * otfp_engine_center_force ( otfp_engine * engine, int center ) {
  if (!engine || !engine->ds) return NULL;
  return DataSpace_centerForce(engine->ds,center);
}

int otfp_engine_add_cv ( otfp_engine * engine, const char * type, int ncenters, const int * centers,
    double zmin, double zmax, const char * boundary, double boundary_k ) {
  if (!engine || !engine->ds || !type || !centers || !boundary) return -1;
  if (!DataSpace_add_cv(engine->ds,(char *)type,ncenters,(int *)centers,
      zmin,zmax,(char *)boundary,boundary_k,"",-1)) return -1;
  return 0;
}

int otfp_engine_add_restraint ( otfp_engine * engine, double k, double target,
    int ncvs, const double * coefficients, const char * function,
    double zmin, double zmax, const char * boundary, double boundary_k ) {
  restraint * r;

  if (!engine || !engine->ds || !coefficients || !function || !boundary) return -1;
  if (engine->irestr>=engine->nrestr) return -1;

  r=DataSpace_AddRestr(engine->ds,k,target,ncvs,(double *)coefficients,
      (char *)function,zmin,zmax,(char *)boundary,boundary_k,"",-1);
  if (!r) return -1;

  engine->restraints[engine->irestr]=r;
  engine->irestr++;
  return engine->irestr-1;
}

static restraint * otfp_engine_restraint ( otfp_engine * engine, int index ) {
  if (!engine || !engine->restraints) return NULL;
  if (index<0 || index>=engine->irestr) return NULL;
  return engine->restraints[index];
}

int otfp_engine_add_tamd ( otfp_engine * engine, int restraint,
    double gamma, double kt, double dt, int chapeau_id, int chapeau_dim ) {
  return restr_AddTamdOpt(otfp_engine_restraint(engine,restraint),
      gamma,kt,dt,chapeau_id,chapeau_dim);
}

int otfp_engine_add_smd ( otfp_engine * engine, int restraint,
    double target, int first_step, int last_step ) {
  return restr_AddSmdOpt(otfp_engine_restraint(engine,restraint),
      target,first_step,last_step);
}

int otfp_engine_add_adam ( otfp_engine * engine, int restraint,
    double step_size, double beta1, double beta2, double epsilon,
    double decay, double kt, double min_step ) {
  return restr_AddAdamOpt(otfp_engine_restraint(engine,restraint),
      step_size,beta1,beta2,epsilon,decay,kt,min_step);
}

int otfp_engine_compute ( otfp_engine * engine, int first_step, int step, double external_bias ) {
  if (!engine || !engine->ds) return -1;
  DataSpace_ComputeCVs(engine->ds);
  return DataSpace_RestrainingForces(engine->ds,first_step,step,external_bias);
}

double otfp_engine_bias_energy ( otfp_engine * engine ) {
  if (!engine || !engine->ds) return 0.0;
  return DataSpace_RestraintEnergy(engine->ds);
}
