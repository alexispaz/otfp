#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "otfp_engine.h"

/* Two Adam restraints that terminate at different steps. The restraint
 * checkpoint must be written when each one terminates, not on every step
 * after the first one ends. */

static const char *ckpt = "checkpoint.rstr";

static double noise(void) {
  double v1, v2, rsq;
  do {
    v1 = 2.0*drand48()-1.0;
    v2 = 2.0*drand48()-1.0;
    rsq = v1*v1+v2*v2;
  } while (rsq >= 1.0 || rsq == 0.0);
  return v1*sqrt(-2.0*log(rsq)/rsq);
}

static int step(otfp_engine *engine, int t) {
  double p0[] = {0.0, 0.0, 0.0};
  double p1[] = {5.0+noise(), 0.0, 0.0};
  otfp_engine_set_center_position(engine, 0, p0);
  otfp_engine_set_center_position(engine, 1, p1);
  return otfp_engine_compute(engine, 0, t, 0.0);
}

int main(void) {
  int ind0[] = {0};
  int ind1[] = {1};
  int cv_centers[] = {0, 1};
  double mass[] = {1.0};
  double cvc[] = {1.0};
  int r, t, code = 0, first_end = 0, writes_after = 0;

  otfp_engine *engine = otfp_engine_create(2, 1, 2, 1234);
  if (!engine) return 1;
  if (otfp_engine_add_center(engine, 1, ind0, mass) != 0) return 2;
  if (otfp_engine_add_center(engine, 1, ind1, mass) != 1) return 3;
  if (otfp_engine_add_cv(engine, "BOND", 2, cv_centers, 0.0, 0.0,
      "NADA", 0.0) != 0) return 4;
  for (r = 0; r < 2; r++) {
    if (otfp_engine_add_restraint(engine, 1.0, 4.0, 1, cvc, "HARMONIC",
        0.0, 0.0, "NADA", 0.0) != r) return 5;
    if (otfp_engine_add_adam(engine, r, 0.01, 0.9, 0.99, 1e-3, 1e-3, 1.0,
        r ? 4000 : 1000) != 0) return 6;
  }

  srand48(5);
  remove(ckpt);
  /* Steps start at 1 and stay below the periodic save frequency, so every
   * checkpoint seen here comes from a termination. */
  for (t = 1; t < 90000 && code != 2; t++) {
    code = step(engine, t);
    if (access(ckpt, F_OK) == 0) {
      if (first_end) writes_after++;
      else first_end = t;
      remove(ckpt);
    }
  }
  otfp_engine_free(engine);

  if (code != 2 || !first_end) {
    fprintf(stderr, "restraints did not both terminate (code %d, first %d)\n",
            code, first_end);
    return 7;
  }
  /* One write when the second restraint terminates, none in between. */
  if (writes_after != 1) {
    fprintf(stderr, "first restraint ended at step %d, all at %d; "
            "%d checkpoint writes in between and at the end, expected 1\n",
            first_end, t-1, writes_after);
    return 8;
  }
  return 0;
}
