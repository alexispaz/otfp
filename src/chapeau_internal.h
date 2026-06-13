#pragma once
#ifndef OTFP_CHAPEAU_INTERNAL_H
#define OTFP_CHAPEAU_INTERNAL_H

#include <stdio.h>

#include "chapeau.h"

#ifdef USE_MKL
#include <mkl.h>
#include <mkl_dss.h>
#endif

// Chapeau type definition: this defines a linear expansion of chapeau
// functions and stores the statistics used to evolve the coefficients.
struct CHAPEAU {
  int nupdate;

  int dm;
  int * N;
  int m;
  int ku;
  int ldad;

  int * periodic;
  double * rmin;
  double * rmax;
  double * dr;
  double * idr;

  double * r;
  double * f;

  FILE * ofp;
  int outputFreq;
  int outputLevel;
  char filename[255];

  double * b;
  double * lam;
  int * hits;

  double norm;

#ifdef USE_MKL
  _MKL_DSS_HANDLE_t * handle;
  MKL_INT nRows, non0;
  int * rowIndex;
  int * columns;
  double * A;
#else
  double ** A;
#endif

  double * bfull;
#ifdef USE_MKL
  double * Afull;
#else
  double ** Afull;
#endif

  int (*accumulate)(struct CHAPEAU * self, double bias);
};

int chapeau_init ( chapeau * ch, int dm, double * rmin, double * rmax, int * N, int * periodic);
int chapeau_comparesize ( chapeau * ch1, chapeau * ch2);
int chapeau_comparegrid ( chapeau * ch1, chapeau * ch2);
void chapeau_set_update_frequency ( chapeau * ch, int nupdate );
void chapeau_set_sample ( chapeau * ch, int dim, double r, double f );
int chapeau_accumulate_sample ( chapeau * ch, double bias );
int chapeau_should_solve ( chapeau * ch, int timestep, int beginsolve );
int chapeau_should_output ( chapeau * ch, int timestep );
void chapeau_savestate_default ( chapeau * ch );

void chapeau_sum ( chapeau * ch1, chapeau * ch2 );
chapeau * chapeau_crop (chapeau * ch, double * rmin, double * rmax);
int chapeau_simetrize2D ( chapeau * ch);

void chapeau_solve ( chapeau * ch );
void chapeau_solve_secure ( chapeau * ch );

double chapeau_evalf_1simplex ( chapeau * ch, double z );
double chapeau_evalf_2simplex ( chapeau * ch, double z1, double z2 );
double chapeau_evalfg_2simplex ( chapeau * ch, double z1, double z2, double * g);
int chapeau_caneval_f1 ( chapeau * ch, double z );
int chapeau_caneval_f2 ( chapeau * ch, double z1, double z2 );
void chapeau_set_peaks ( chapeau * ch, char * filename );

char * chapeau_serialize ( chapeau * ch );
void chapeau_addserialized ( chapeau * ch, char * str );
void chapeau_setserialized ( chapeau *ch, char * str );

int accumulate_1D ( chapeau * ch, double bias);
int accumulate_2D ( chapeau * ch, double bias );

#endif
