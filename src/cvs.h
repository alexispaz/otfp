#pragma once
#ifndef OTFP_CVS_H
#define OTFP_CVS_H

typedef struct cv cv;

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
