#include <stdio.h>

#include "otfp_engine.h"

int main(void) {
  int atom0[] = {0};
  int atom1[] = {1};
  int cv_centers[] = {0, 1};
  double mass[] = {1.0};
  double coefficient[] = {1.0};
  double p0[] = {0.0, 0.0, 0.0};
  double p1[] = {1.0, 0.0, 0.0};
  const double *f0;
  const double *f1;
  int restraint;
  int status;
  int exit_code = 0;

  otfp_engine *engine = otfp_engine_create(2, 1, 1, 1234);
  if (!engine) return 1;

  if (otfp_engine_add_center(engine, 1, atom0, mass) != 0) { exit_code = 2; goto done; }
  if (otfp_engine_add_center(engine, 1, atom1, mass) != 1) { exit_code = 3; goto done; }
  if (otfp_engine_add_cv(engine, "BOND", 2, cv_centers, 0.0, 0.0,
      "NADA", 0.0) != 0) { exit_code = 4; goto done; }

  restraint = otfp_engine_add_restraint(engine, 1.0, 0.0, 1, coefficient,
      "HARMONIC", 0.0, 0.0, "NADA", 0.0);
  if (restraint < 0) { exit_code = 5; goto done; }
  if (otfp_engine_add_smd(engine, restraint, 1.0, 0, 10) != 0) { exit_code = 6; goto done; }

  if (otfp_engine_set_center_position(engine, 0, p0) != 0) { exit_code = 7; goto done; }
  if (otfp_engine_set_center_position(engine, 1, p1) != 0) { exit_code = 8; goto done; }

  status = otfp_engine_compute(engine, 0, 1, 0.0);
  if (status != 0) { exit_code = 9; goto done; }

  f0 = otfp_engine_center_force(engine, 0);
  f1 = otfp_engine_center_force(engine, 1);

  printf("energy %.12f\n", otfp_engine_bias_energy(engine));
  printf("force 0 %.12f %.12f %.12f\n", f0[0], f0[1], f0[2]);
  printf("force 1 %.12f %.12f %.12f\n", f1[0], f1[1], f1[2]);

done:
  otfp_engine_free(engine);
  return exit_code;
}
