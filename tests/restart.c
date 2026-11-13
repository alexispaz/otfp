#include <math.h>
#include <stdio.h>

#include "cfacv.h"

static int setup_dataspace(DataSpace *ds, double target, restraint **out) {
  int ind0[] = {0};
  int ind1[] = {1};
  int cv_ind[] = {0, 1};
  double mass[] = {1.0};
  double cvc[] = {1.0};

  if (DataSpace_AddAtomCenter(ds, 1, ind0, mass) != 0) return 1;
  if (DataSpace_AddAtomCenter(ds, 1, ind1, mass) != 1) return 2;
  if (!DataSpace_add_cv(ds, "BOND", 2, cv_ind, 0.0, 0.0, "NADA", 0.0, "", -1)) return 3;

  *out=DataSpace_AddRestr(ds, 1.0, target, 1, cvc, "HARMONIC",
      0.0, 0.0, "NADA", 0.0, "", -1);
  if (!*out) return 4;

  return 0;
}

int main(void) {
  const char *filename = "restart-test.rstr";
  DataSpace *src = NewDataSpace(2, 1, 1, 1234);
  DataSpace *dst = NewDataSpace(2, 1, 1, 1234);
  restraint *src_r = NULL;
  restraint *dst_r = NULL;
  int code;

  if (!src || !dst) return 1;

  code=setup_dataspace(src, 2.5, &src_r);
  if (code) return 10+code;

  code=setup_dataspace(dst, 0.0, &dst_r);
  if (code) return 20+code;

  DataSpace_save_restraints(src, (char *)filename);
  DataSpace_load_restraints(dst, (char *)filename);

  remove(filename);

  if (fabs(restr_getz(src_r)-restr_getz(dst_r))>1e-12) return 30;
  if (fabs(restr_getz(dst_r)-2.5)>1e-12) return 31;

  /* The first step of the restarted run must keep the loaded target, not
   * reset it to the CV value (1.0 here). */
  {
    double *x0=DataSpace_centerPosition(dst,0);
    double *x1=DataSpace_centerPosition(dst,1);
    x0[0]=0.0; x0[1]=0.0; x0[2]=0.0;
    x1[0]=1.0; x1[1]=0.0; x1[2]=0.0;
    DataSpace_ComputeCVs(dst);
    DataSpace_RestrainingForces(dst,1,1,0.0);
    if (fabs(restr_getz(dst_r)-2.5)>1e-12) {
      fprintf(stderr,"first step reset the loaded target to %g\n",restr_getz(dst_r));
      return 32;
    }
  }

  /* Without a restart, the first step still starts the target at the CV. */
  {
    double *x0=DataSpace_centerPosition(src,0);
    double *x1=DataSpace_centerPosition(src,1);
    x0[0]=0.0; x0[1]=0.0; x0[2]=0.0;
    x1[0]=1.0; x1[1]=0.0; x1[2]=0.0;
    DataSpace_ComputeCVs(src);
    DataSpace_RestrainingForces(src,1,1,0.0);
    if (fabs(restr_getz(src_r)-1.0)>1e-12) return 33;
  }

  DataSpace_free(src);
  DataSpace_free(dst);

  return 0;
}
