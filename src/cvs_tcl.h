#pragma once
#ifndef OTFP_CVS_TCL_H
#define OTFP_CVS_TCL_H

typedef struct CHAPEAU chapeau;
typedef struct cv cv;

chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int * periodic );
void chapeau_free ( chapeau * ch );
void chapeau_setupoutput ( chapeau * ch, char * outfile, char * restartfile, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );
void chapeau_savestate ( chapeau * ch, char * filename );
chapeau * chapeau_allocloadstate ( char * filename );
void chapeau_loadstate ( chapeau * ch, char * filename );
void chapeau_loadlambda ( chapeau * ch, char * filename );

cv * cv_init ( char * typ, int nC, int * ind,
    double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);
void cv_free ( cv * c );

int set_zsd_circle ( double x, double y, double xy, double s );
int set_zsd_ring ( double x, double y, double r1, double r2, double s );
int set_d2chap ( cv * c, char * filename );
int set_line ( cv * c );

double * cv_access_ref ( cv * c, int i );
double * cv_access_ref2 ( cv * c, int i );

void cv_output ( cv * c );

#endif
