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

struct DataSpace {
  int N,iN;
  int K,iK;
  int ch_num,ch_now;

  atomCenterStruct ** ac;
  restraint ** restr;
  chapeau ** ch;

  int ncv,icv;
  cv ** cv;

  char filename[16];
  int restrsavefreq;

  double ** R;

  int doAnalyticalCalc;

  double lamfric;
  double lamdt;

  int evolveAnalyticalParameters;
  int beginaccum;
  int beginsolve;
};

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

void ds_saverestrains ( DataSpace * ds, char * filename );
void ds_loadrestrains ( DataSpace * ds, char * filename );

#endif
