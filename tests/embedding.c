#include <math.h>
#include <stdio.h>

#include "otfp_engine.h"

static int is_finite3(const double *v) {
  return v && isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
}

int main(void) {
  int ind0[] = {0};
  int ind1[] = {1};
  int cv_centers[] = {0, 1};
  double mass[] = {1.0};
  double cvc[] = {1.0};
  double p0[] = {0.0, 0.0, 0.0};
  double p1[] = {1.0, 0.0, 0.0};
  const double *f0;
  const double *f1;
  double energy;
  int code;

  otfp_engine *engine = otfp_engine_create(2, 1, 1, 1234);
  if (!engine) {
    fprintf(stderr, "otfp_engine_create returned NULL\n");
    return 1;
  }

  if (otfp_engine_add_center(engine, 1, ind0, mass) != 0) return 2;
  if (otfp_engine_add_center(engine, 1, ind1, mass) != 1) return 3;
  if (otfp_engine_add_cv(engine, "BOND", 2, cv_centers, 0.0, 0.0,
      "NADA", 0.0) != 0) return 4;
  if (otfp_engine_add_restraint(engine, 1.0, 0.0, 1, cvc, "HARMONIC",
      0.0, 0.0, "NADA", 0.0) != 0) return 5;

  if (otfp_engine_set_center_position(engine, 0, p0) != 0) return 6;
  if (otfp_engine_set_center_position(engine, 1, p1) != 0) return 7;

  code = otfp_engine_compute(engine, 0, 1, 0.0);
  if (code != 0) return 8;

  f0 = otfp_engine_center_force(engine, 0);
  f1 = otfp_engine_center_force(engine, 1);
  energy = otfp_engine_bias_energy(engine);

  if (!isfinite(energy)) return 9;
  if (!is_finite3(f0)) return 10;
  if (!is_finite3(f1)) return 11;

  otfp_engine_free(engine);
  otfp_engine_free(NULL);

  return 0;
}
