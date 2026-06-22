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

  DataSpace_free(src);
  DataSpace_free(dst);

  return 0;
}
