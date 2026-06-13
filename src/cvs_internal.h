#pragma once
#ifndef OTFP_CVS_INTERNAL_H
#define OTFP_CVS_INTERNAL_H

#include <stdio.h>

#include "chapeau.h"
#include "cvs.h"

struct cv {
  int typ;
  int nC;
  int * ind;
  double val;
  double ** gr;

  double f;
  double u;
  int amd;

  double min;
  double max;
  double half_domain;
  double boundk;
  int (*boundFunc)(struct cv * self);

  double ** ref;
  double ** ref2;
  double refmod;

  FILE * ofp;
  int outputFreq;
  int boutput;

  int (*calc)(struct cv * self, double ** R );
  chapeau * d2ch;
};

int calccv_line  ( cv * c, double ** R );
int calccv_rmsd  ( cv * c, double ** R );
int calccv_charmmd ( cv * c, double ** R );
int calccv_d2chap  ( cv * c, double ** R );
int set_d2chap    ( cv * c, char * filename );
int calccv_gauss2d ( cv * c, double ** R );
int calccv_cogx  ( cv * c, double ** R );
int calccv_cogy  ( cv * c, double ** R );
int calccv_cogz  ( cv * c, double ** R );
int calccv_x     ( cv * c, double ** R );
int calccv_y     ( cv * c, double ** R );
int calccv_z     ( cv * c, double ** R );
int calccv_s     ( cv * c, double ** R );
int calccv_halfbond  ( cv * c, double ** R );
int calccv_bond  ( cv * c, double ** R );
int calccv_bonds ( cv * c, double ** R );
int calccv_dihed ( cv * c, double ** R );
int calccv_angle ( cv * c, double ** R );
int calccv_zsd_circle ( cv * c, double ** R );
int calccv_zsd_xrange ( cv * c, double ** R );
int calccv_zsd_ring ( cv * c, double ** R );

double cdf(double x);

int cv_amd ( cv * c );
int cv_SoftUpperWall ( cv * c );
int cv_SoftLowerWall ( cv * c );
int cv_SoftWalls ( cv * c );
int cv_nada ( cv * c );

int set_line (cv * c);

int cv_compute ( cv * c, double ** R );
double cv_bias_value ( cv * c );
double cv_value ( cv * c );
int cv_center_count ( cv * c );
int cv_center_index ( cv * c, int i );
double cv_gradient_component ( cv * c, int center, int dim );
void cv_reset_force_energy ( cv * c );
void cv_apply_boundary ( cv * c );
double cv_force ( cv * c );
int cv_should_output ( cv * c, int timestep );

#endif
