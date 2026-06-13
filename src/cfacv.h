#pragma once
#ifndef OTFP_CFACV_H
#define OTFP_CFACV_H

#include <stdio.h>

#include "chapeau.h"
#include "cvs.h"

typedef struct TamdOptStruct tamdOptStruct;
typedef struct SMDOPTSTRUCT smdOptStruct;
typedef struct AdamOptStruct adamOptStruct;
typedef struct restraint restraint;
typedef struct DataSpace DataSpace;

void cfacvBanner ( void );

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp);
smdOptStruct * New_smdOptStruct ( double target, int t0, int t1);
adamOptStruct * New_adamOptStruct ( double a, double b1, double b2, double e , double decay, double kkT, double mvar, double vvar, double tmin);
int smdOptInit ( smdOptStruct * smd, double initval);
int rf_getityp ( char * typ );
char * rf_getstyp ( int ityp );

restraint * New_restraint ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax, char * boundstr, double boundk, char * outfile, int outputFreq);
void restraint_free ( restraint * r );

FILE * my_fopen ( char * name, char * code ) ;
DataSpace * NewDataSpace ( int Nctr, int Ncv, int Nrstr, long int seed );
void DataSpace_free ( DataSpace * ds );
int DataSpace_SetupChapeau ( DataSpace * ds, int numrep, int dm, double * min,
    int * nKnots, double * max, int * periodic, int beginaccum, int beginsolve,
    char * outfile, int outfreq, int outlevel, int nupdate);
chapeau * DataSpace_get_chapeauadress ( DataSpace * ds, int i );
int DataSpace_getN ( DataSpace * ds );
double * DataSpace_centerPos ( DataSpace * ds, int i );
int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
cv * DataSpace_add_cv ( DataSpace * ds, char * typ, int nind, int * ind,
    double zmin, double zmax,char * boundf, double boundk, char * outfile, int outputFreq );

restraint * DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax,char * boundf, double boundk,char * outfile, int outputFreq);
int restr_UpdateTamdOpt ( restraint * r, double g, double kt, double dt );
int restr_AddTamdOpt ( restraint * r, double g, double kt, double dt, int chid, int chdm );
int restr_AddSmdOpt  ( restraint * r, double target, int t0, int t1 );
int restr_AddAdamOpt  ( restraint * r, double a, double b1, double b2, double e, double decay, double kT, double tmin);
void restr_output  ( restraint * r );
double DataSpace_ComputeCVs ( DataSpace * ds);
int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep, double bias);
double DataSpace_RestraintEnergy ( DataSpace * ds );
void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
int DataSpace_checkdata ( DataSpace * ds );
int DataSpace_dump ( DataSpace * ds );
FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );

double restr_getz ( restraint * r );
double restr_getu ( restraint * r );
int restr_set_rchid ( restraint * r, DataSpace * ds, int chid);

#endif
