#include <stdio.h>

#include "cfacv.h"
#include "chapeau.h"
#include "cvs.h"

int main(void) {
  int ind0[] = {0};
  int ind1[] = {1};
  int cv_ind[] = {0, 1};
  double mass[] = {1.0};
  double cvc[] = {1.0};

  DataSpace *ds = NewDataSpace(2, 1, 1, 1234);
  if (!ds) {
    fprintf(stderr, "NewDataSpace returned NULL\n");
    return 1;
  }

  if (DataSpace_AddAtomCenter(ds, 1, ind0, mass) != 0) return 2;
  if (DataSpace_AddAtomCenter(ds, 1, ind1, mass) != 1) return 3;

  if (!DataSpace_add_cv(ds, "BOND", 2, cv_ind, 0.0, 0.0, "NADA", 0.0, "", -1)) {
    return 4;
  }

  if (!DataSpace_AddRestr(ds, 1.0, 0.0, 1, cvc, "HARMONIC", 0.0, 0.0, "NADA", 0.0, "", -1)) {
    return 5;
  }

  DataSpace_free(ds);
  DataSpace_free(NULL);
  restraint_free(NULL);
  cv_free(NULL);
  chapeau_free(NULL);

  return 0;
}
