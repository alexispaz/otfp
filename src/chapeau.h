#pragma once
#ifndef OTFP_CHAPEAU_H
#define OTFP_CHAPEAU_H

typedef struct CHAPEAU chapeau;

chapeau * chapeau_alloc ( int dm, double * rmin, double * rmax, int * N, int * periodic );
void chapeau_free ( chapeau * ch );

void chapeau_setupoutput ( chapeau * ch, char * outfile, char * restartfile, int outputFreq, int outputLevel );
void chapeau_output ( chapeau * ch, int timestep );

void chapeau_savestate ( chapeau * ch, char * filename );
chapeau * chapeau_allocloadstate ( char * filename );
void chapeau_loadstate ( chapeau * ch, char * filename );
void chapeau_loadlambda ( chapeau * ch, char * filename );

#endif
