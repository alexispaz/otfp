#pragma once
#ifndef OTFP_CFACV_INTERNAL_H
#define OTFP_CFACV_INTERNAL_H

#include <stdio.h>

#include "centers.h"
#include "cfacv.h"
#include "chapeau.h"
#include "cvs.h"

struct TamdOptStruct {
  double kT;
  double gamma;
  double dt;
  double ginv;
  double noise;
};

typedef double (*smdUpdateFunc) (double, double, int );
struct SMDOPTSTRUCT {
  int t0;
  int t1;
  double invinterval;
  double initval;
  double increment;
  double target;
  smdUpdateFunc update;
};

struct AdamOptStruct {
  int t;
  int tmin;
  int stage;
  double e;
  double b1;
  double b2;
  double a;
  double decay;
  double v;
  double m;
  double mvar;
  double vvar;
  double vmed;
  double zvar;
};

enum {HARMONIC, HARMCUTO, NULL_RF};

struct restraint {
  double k;
  double z;
  double val;

  double f;
  double u;
  int nCV;
  double * cvc;

  int evolve;
  int restored;

  tamdOptStruct * tamdOpt;
  double tamd_noise;
  double tamd_restraint;

  smdOptStruct * smdOpt;
  adamOptStruct * adamOpt;

  int rfityp;
  int (*evolveFunc)(struct restraint * self);
  int (*energyFunc)(struct restraint * self);

  double min;
  double max;
  double half_domain;
  double boundk;
  int (*boundFunc)(struct restraint * self);

  int chid;
  int chdm;

  FILE * ofp;
  int outputFreq;
  int boutput;
};

typedef struct CenterStepState {
  int ctr_count;
  double ** ctr_positions;
  double ** ctr_forces;
} CenterStepState;

typedef struct CenterSet {
  int ctr_total;
  int ctr_count;
  atomCenterStruct ** ctrs;
} CenterSet;

typedef struct CvSet {
  int total;
  int count;
  cv ** cvs;
} CvSet;

typedef struct RestraintSet {
  int total;
  int count;
  restraint ** items;
  char checkpoint_file[16];
  int save_frequency;
} RestraintSet;

typedef struct FesState {
  int total;
  int enabled;
  int evolve_parameters;
  int begin_accum;
  int begin_solve;
  chapeau ** chapeaux;
} FesState;

struct DataSpace {
  CenterSet centers;
  CvSet cvs;
  RestraintSet restraints;
  CenterStepState center_step;
  FesState fes;
};

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp);
smdOptStruct * New_smdOptStruct ( double target, int t0, int t1);
adamOptStruct * New_adamOptStruct ( double a, double b1, double b2, double e , double decay, double kkT, double mvar, double vvar, double tmin);
int smdOptInit ( smdOptStruct * smd, double initval);
int rf_getityp ( char * typ );
char * rf_getstyp ( int ityp );

restraint * New_restraint ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);
void restraint_free ( restraint * r );
void restr_output  ( restraint * r );

FILE * my_fopen ( char * name, char * code ) ;
FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );

int DataSpace_getN ( DataSpace * ds );
int DataSpace_checkdata ( DataSpace * ds );
int DataSpace_dump ( DataSpace * ds );
int DataSpace_ctr_count ( DataSpace * ds );
int DataSpace_cv_count ( DataSpace * ds );
int DataSpace_rstr_count ( DataSpace * ds );
int DataSpace_fes_count ( DataSpace * ds );
int DataSpace_fes_is_enabled ( DataSpace * ds );
int DataSpace_fes_evolves ( DataSpace * ds );
int DataSpace_fes_accum_begin ( DataSpace * ds );
int DataSpace_fes_solve_begin ( DataSpace * ds );
int DataSpace_rstr_save_frequency ( DataSpace * ds );
char * DataSpace_rstr_checkpoint_file ( DataSpace * ds );
cv * DataSpace_cv ( DataSpace * ds, int index );
restraint * DataSpace_rstr ( DataSpace * ds, int index );
chapeau * DataSpace_fes ( DataSpace * ds, int index );
double ** DataSpace_ctr_positions ( DataSpace * ds );
double * DataSpace_ctr_position ( DataSpace * ds, int ctr );
double * DataSpace_ctr_force ( DataSpace * ds, int ctr );
void DataSpace_clear_ctr_forces ( DataSpace * ds );
void DataSpace_add_ctr_force ( DataSpace * ds, int ctr, int dim, double force );

int cbd ( restraint * r );
int uniformvelocity ( restraint * r );
int adam ( restraint * r );

int SoftUpperWall ( restraint * r );
int SoftLowerWall ( restraint * r );
int SoftWalls ( restraint * r );
int pbc ( restraint * r );
int Periodic ( restraint * r );
int nada ( restraint * r );

int HarmonicCart ( restraint * r );
int HarmonicCart_cutoff ( restraint * r );
int HarmonicCart_pbc ( restraint * r );
int HarmonicCart_cutoff_pbc ( restraint * r );

void restraint_reset_value ( restraint * r );
void restraint_add_cv_value ( restraint * r, int cv_index, double value );
void restraint_init_target_from_value ( restraint * r );
int restraint_compute_energy ( restraint * r );
int restraint_should_output ( restraint * r, int timestep );
int restraint_cv_count ( restraint * r );
double restraint_cv_coefficient ( restraint * r, int cv_index );
double restraint_force ( restraint * r );
double restraint_target ( restraint * r );
int restraint_chapeau_id ( restraint * r );
int restraint_chapeau_dim ( restraint * r );
int restraint_evolve_for_timestep ( restraint * r, int timestep );
void restraint_stop_evolution ( restraint * r );
int restraint_is_evolving ( restraint * r );
void restraint_report ( restraint * r, int index );

void DataSpace_save_restraints ( DataSpace * ds, char * filename );
void DataSpace_load_restraints ( DataSpace * ds, char * filename );

#endif
