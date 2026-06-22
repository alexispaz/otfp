/* Collective Variables module for NAMD via tclForces
 * Cameron F. Abrams
 *
 * Extended for 2D OTFP and RE OTFP
 * Sergio A. Paz
 *
 * This header declares structured data types and functions for
 * handling collective variables and their restraints.
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
 * transmit them to the atoms. This two-way communication happens
 * once per timestep on the master node during NAMD execution.
 *
 * In this particular OTFP impelementation, we are using linear chapeau
 * functions as the basis functions. We use this to compute effective pair
 * potentials between centers of mass of atomgroups (like water molecules).
 *
 */

#include <stdio.h>
#include <string.h>

#include "cfacv_internal.h"
#include "cvs_internal.h"
#include "cfacv.h"

void cfacvBanner ( void ) {
  printf("OTFP: Version (Git Hash): %s\n",VERSION);
}

FILE * my_fopen ( char * name, char * code ) {
  if (strcmp(name,"stdout")) return fopen(name,code);
  else return stdout;
}

FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds ) {
  FILE * fp=fopen(name,code);
  fwrite(&outputLevel,sizeof(unsigned int),1,fp);
  {
    int nrestr=DataSpace_rstr_count(ds);
    fwrite(&nrestr,sizeof(int),1,fp);
  }
  return fp;
}


double DataSpace_ComputeCVs ( DataSpace * ds ) {
  int i;
  double bias;
 
  bias=0;
  if (!ds) return -1;
  for (i=0;i<DataSpace_cv_count(ds);i++){
    cv * c=DataSpace_cv(ds,i);
    cv_compute(c,DataSpace_ctr_positions(ds));
    bias+=cv_bias_value(c);
  }
 
  return bias;
}

double DataSpace_RestraintEnergy ( DataSpace * ds ) {
  if (ds) {
    double u=0.0;
    int i;
    restraint * r;
    for (i=0;i<DataSpace_rstr_count(ds);i++) {
      r=DataSpace_rstr(ds,i);
      u+=restr_getu(r);
    }
    return u;
  }
  return 0.0;    
}

void DataSpace_ReportCV ( DataSpace * ds, int * active, double * res ) {
  if (ds) {
    int i;
    for (i=0;i<DataSpace_cv_count(ds);i++) {
      if (active[i]) res[i]=cv_value(DataSpace_cv(ds,i));
      else res[i]=0.0;
    }
  }
}

void DataSpace_ReportAll ( DataSpace * ds ) {
  if (ds) {
    int i;
    printf("CFACV_OTFP/C) DataSpace: %i Centers, %i CV's, %i Restraints\n",
	   DataSpace_ctr_count(ds), DataSpace_cv_count(ds), DataSpace_rstr_count(ds));
    for (i=0;i<DataSpace_ctr_count(ds);i++) {
      double * ctr = DataSpace_ctr_position(ds,i);
      printf("CFACV_OTFP/C)  Center %i : %.5f %.5f %.5f\n",
	     i,ctr[0],ctr[1],ctr[2]);
    }
    for (i=0;i<DataSpace_cv_count(ds);i++) {
      cv_report(DataSpace_cv(ds,i),i);
    }
    for (i=0;i<DataSpace_rstr_count(ds);i++) {
      restraint_report(DataSpace_rstr(ds,i),i);
    }
    
  }
}
