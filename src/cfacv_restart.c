#include <stdio.h>
#include <stdlib.h>

#include "cfacv_internal.h"

void DataSpace_save_restraints ( DataSpace * ds, char * filename ) {
  int i;
  restraint * r;
  FILE * ofs;

  ofs=fopen(filename,"w");
  {
    int nrestr=DataSpace_rstr_count(ds);
    fwrite(&nrestr,sizeof(int),1,ofs);
  }
  for (i=0;i<DataSpace_rstr_count(ds);i++) {
    r=DataSpace_rstr(ds,i);
    fwrite(&(r->z),sizeof(double),1,ofs);
    fwrite(&(r->val),sizeof(double),1,ofs);
  }

  fclose(ofs);
}

void DataSpace_load_restraints ( DataSpace * ds, char * filename ) {
  int i;
  restraint * r;
  FILE * ofs;

  ofs=fopen(filename,"r");

  fread(&i,sizeof(int),1,ofs);
  fprintf(stdout,"CFACV) Recovering restraint state of step %i of some previous simulation\n",i);
  if ( i!=DataSpace_rstr_count(ds) ) {
    fprintf(stderr,"CFACV) ERROR restraint restart because holding object has not the proper size\n");
    exit(-1);
  }

  for (i=0;i<DataSpace_rstr_count(ds);i++) {
    r=DataSpace_rstr(ds,i);
    fread(&(r->z),sizeof(double),1,ofs);
    fread(&(r->val),sizeof(double),1,ofs);
    printf("OTFP: Recovering restraint %i, z=%.5f; val=%.5f \n",i,r->z,r->val);
  }

  fclose(ofs);
}
