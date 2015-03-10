#pragma once

/* Collective Variables module for NAMD via tclForces
 * Cameron F. Abrams
 *
 * Adapted for On-the-fly parameterization of free energies via TAMD
 *
 * This header declares structured data types and functions for
 * handling collective variables and their restraints.  These
 * functions are defined in cfacv.c and invoked in cfacv.tcl.
 * cfacv.tcl is the master tcl script that defines all the tcl
 * functions used by the tclforces interface script,
 * cfacv_tclforces.tcl.  cfacv_tclforces.tcl must be referenced in
 * the NAMD configuration file with the "tclforcesscript" keyword.
 * 
 * The central concept is that of a "restraint" which is applied to
 * one or a linear combination of "collective variables".  A restraint
 * computes a force on a collective variable ($\theta$) based on its
 * value and the value of a "target" or "auxiliary" variable ($z$).
 * The auxiliary variable can be constant, or driven according to
 * steered MD or temperature-accelerated MD.  It is the job of the
 * cfacv module to compute these forces after receiving from the NAMD
 * tclforces "calcforces" function the positions of all atoms or
 * atom-groups that participate in each CV owned by each restraint.
 * It is the job of the tclforces module to receive those forces and
 * transmit them to the atoms.  This two-way communication happens
 * once per timestep on the master node during NAMD execution.
 *
 * The most important data type defined here is the "dataspace", which
 * is just the working space that receives the atom/atom-group
 * positions, computes the values of the CV's, the restraints, the
 * restraint auxiliaries, and the resultant forces on each atom or
 * atom-group.
 *
 * In this version of OTFP, the free energy is a sum of effective
 * pair potentials:
 *
 * G(z;\lambda) = \sum_{i<j} g(|z_i-z_j|;\lambda)
 *
 * where
 * 
 * g(z;\lambda) = \sum_k \lambda_k \phi_k(z)
 *
 * is a pair potential cast as a linear expansion in basis functions
 * \phi_k(z).  Generally, z is the set of variables auxiliary to the
 * monitored CV's, and here, these correspond to the cartesian
 * positions of each atom or atom-block considered.  In this particular
 * OTFP impelementation, we are using linear chapeau functions as 
 * the basis functions.  We use this to compute effective pair potentials
 * between centers of mass of atomgroups (like water molecules).
 * 
 * Data in Fig. 5 of this paper were computed using this code: 
 *
 * Abrams and Vanden-Eijnden, "On-the-fly free energy parameterization
 * via temperature-accelerated molecular dynamics," Chem Phys Lett
 * 547:114-119 (2012).
 *
 * Modifications to this code should allow computation of 1-D free
 * energies for an arbitrary collective variable.  Slightly more
 * involved modifications should allow for computation of
 * multidimensional free energies as multilinear basis function
 * expansions of the form
 *
 * G(z;\lambda) = \sum_{klm...} \lambda_{klm...}\phi_k(z_1)\phi_l(z_2)\phi_m(z_3)...
 *
 * c 2009-2014 Cameron F Abrams
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "centers.h"
//#include "wrapcoords.h"
#include "chapeau.h"       // basis functions

// Cross reference with cfacv.h
struct cvStruct;  
typedef struct CVSTRUCT cvStruct;
#include "cvs.h"           // collective variables


// Temperature-accelerated MD options structure; will be owned by the
// restraint structure
typedef struct TAMDOPTSTRUCT {
  double kT;           // fictitious thermal energy (kcal/mol)
  double gamma;        // fictitious friction (1/ps)
  double dt;           // timestep (ps)
  double ginv;         // inverse friction
  double noise;        // noise
  int periodic;        // flag indicating whether restraint is
                       // periodic (like a dihedral)
  double half_domain;  // half-domain size of periodic CV
} tamdOptStruct;


// Steered MD options structure; will be owned by the restraint
// structure
typedef double (*smdUpdateFunc) (double, double, int );
typedef struct SMDOPTSTRUCT {
  int t0;              // timestep to begin steering
  int t1;              // timestep to end steering (and hold restraint fixed)
  double invinterval;  // 1/(t1-t0)
  double initval;      // initial value of steered variable
  double increment;    // increment of steered variable per timestep
  double target;       // target value of steered variable
  smdUpdateFunc update; // function used to update steered variable
} smdOptStruct;


// The Restraint structure
enum {HARMONIC, HARMCUTO, PERIODIC, NULL_RF};

typedef struct RESTRSTRUCT {
  double k;               // spring constant
  double z;               // target value
  double val;             // restraint value
  double min;             // min, max, half_domain for periodic boundaries
  double max;
  double half_domain;
  double f;              // force on restraint
  double u;              // potential energy stored by restraint
  int nCV;               // number of collective variables in SYSTEM
  double * cvc;          // the coefficient each collective variable
                         // contributes to the restraint.  This allows
                         // a restraint to be applied on any linear
                         // combination of collective variables.  Most
                         // commonly, however, this will be an array
                         // with only one non-zero value (of 1) that
                         // indicates which CV this restraint is
                         // applied to.

  tamdOptStruct * tamdOpt; // pointer to the tamd options structure
  double tamd_noise;
  double tamd_restraint;
  smdOptStruct * smdOpt;   // pointer to the smd options structure

  int rfityp;                   // type of the restraining function
                                // (Harmonic or Periodic)

  int tamd_evolve;        // Indicate when z must evolve. When tamd_evolve goes
                          // from 0 to 1 z value is initialized

  // pointer to the restraining energy and force function
  //restrEnergyFunc energyFunc;   
  int (*energyFunc)(struct RESTRSTRUCT * self, int pbc, double half_domain);

} restrStruct;


void cfacvBanner ( void );

cvStruct * New_cvStruct ( int typ, int nC, int * ind );

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp, double half_domain );

smdOptStruct * New_smdOptStruct ( double target, int t0, int t1, int periodic );

restrStruct * New_restrStruct ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax  );

#ifndef MAXNBOR
#define MAXNBOR 250
#endif


typedef struct DATASPACESTRUCT {
  int N; // number of centers
  int M; // number of CVs
  int K; // number of restraints
  int iN;
  int iM;
  int iK;

  atomCenterStruct ** ac; // defined in centers.h
  cvStruct ** cv;
  restrStruct ** restr;
  chapeau ** ch;

  double ** R; // array of center cartesian coordinates R[i][0/1/2]

  // below are bits for doing TAMD/OTFP of pairwise intercenter
  // potentials represented as piecewise-continuous linear functions

  // PBC
  int pbc;
  double O[3]; // simulation box origin
  double L[3]; // simulation box size
  double hL[3]; // simulation box half-size
  double Min[3], Max[3];

  double squaredPairCutoff;
  int doAnalyticalCalc; // indicates whether or not we are to perform
			// the optmization of the analytical
			// parameterization
                        
  // De alguna manera cada chapeau es una direccion
  int ch_num;
  int * ch_id;

  // for evolution of the parameterization
  double lamfric;
  double lamdt;
  int nsamples;
  int evolveAnalyticalParameters;
  int beginEvolveParameters;
  int useTAMDforces; // indicates whether tamd forces are used to
		     // update auxiliary variables; other choice is
		     // the use forces available from the
		     // instantaneous parameterization of the
		     // analytical free energy (this would make the
		     // calculation conform to the heterogeneous
		     // multiscale method
  int reportParamFreq;
  
} DataSpace;



//The random seed
unsigned short * Xi;

// Other subroutines
FILE * my_fopen ( char * name, char * code ) ;
int * DataSpace_chid ( DataSpace * ds );
DataSpace * NewDataSpace ( int N, int M, int K, long int seed );
int DataSpace_SetupPBC ( DataSpace * ds, int pbc, double Ox, double Oy, double Oz, double Lx, double Ly, double Lz );
int DataSpace_SetupPairCalc ( DataSpace * ds, double cutoff, double nlcutoff, int beginEvolve, int useTAMDforces, int reportParamFreq, double spline_min, int nKnots, char * splineoutputfile, int splineoutputfreq, int splineoutputlevel, int lamupdateinterval, int chnum );
int DataSpace_getN ( DataSpace * ds );
double * DataSpace_centerPos ( DataSpace * ds, int i );
int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) ;
int DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax );
int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt );
int DataSpace_AddSmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 );
int DataSpace_ComputeCVs ( DataSpace * ds );
int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
double DataSpace_RestraintEnergy ( DataSpace * ds );
void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res );
void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
int DataSpace_SetRestraints ( DataSpace * ds, double * rval );
int DataSpace_checkdata ( DataSpace * ds );
int DataSpace_dump ( DataSpace * ds ); 
FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
int DataSpace_InitKnots ( DataSpace * ds, char * filename, int j);



int fes1D( DataSpace * ds ); 
int fes_from_distances( DataSpace * ds, int first, int timestep ) ; 


