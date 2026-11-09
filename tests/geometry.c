#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "measurements.h"

/* Dihedrals and angles at their edges (0 and pi) must stay finite.
 * Rounding can push the cosine just past +-1, and acos() of that is NaN;
 * in a NAMD run the NaN reaches addforce as "force not a vector". */

static void rotation(double R[3][3]) {
  double q[4], n = 0.0, a, b, c, d;
  int i;
  for (i = 0; i < 4; i++) { q[i] = 2.0*drand48()-1.0; n += q[i]*q[i]; }
  n = sqrt(n);
  a = q[0]/n; b = q[1]/n; c = q[2]/n; d = q[3]/n;
  R[0][0] = a*a+b*b-c*c-d*d; R[0][1] = 2*(b*c-a*d);     R[0][2] = 2*(b*d+a*c);
  R[1][0] = 2*(b*c+a*d);     R[1][1] = a*a-b*b+c*c-d*d; R[1][2] = 2*(c*d-a*b);
  R[2][0] = 2*(b*d-a*c);     R[2][1] = 2*(c*d+a*b);     R[2][2] = a*a-b*b-c*c+d*d;
}

static void place(double loc[][3], double P[][3], int n) {
  double R[3][3], o[3];
  int a, i, j;
  rotation(R);
  for (i = 0; i < 3; i++) o[i] = 20.0*(drand48()-0.5);
  for (a = 0; a < n; a++)
    for (i = 0; i < 3; i++) {
      P[a][i] = o[i];
      for (j = 0; j < 3; j++) P[a][i] += R[i][j]*loc[a][j];
    }
}

static int finite3(double g[3]) {
  return isfinite(g[0]) && isfinite(g[1]) && isfinite(g[2]);
}

int main(void) {
  const double edges[] = {0.0, M_PI};
  const long n = 200000;
  long t, bad_dihed = 0, bad_angle = 0;
  int e;

  srand48(11);
  for (e = 0; e < 2; e++) {
    for (t = 0; t < n; t++) {
      double phi = edges[e]+(drand48()-0.5)*2e-8;
      double loc[4][3] = {{-0.5, 1.4, 0.0}, {0.0, 0.0, 0.0}, {1.53, 0.0, 0.0},
                          {2.03, 1.4*cos(phi), 1.4*sin(phi)}};
      double P[4][3], g[4][3], val;
      place(loc, P, 4);
      val = my_getdihed(P[0], P[1], P[2], P[3], g[0], g[1], g[2], g[3]);
      if (!isfinite(val) || fabs(fabs(val)-edges[e]) > 1e-6
          || !finite3(g[0]) || !finite3(g[1]) || !finite3(g[2]) || !finite3(g[3]))
        bad_dihed++;
    }
  }

  /* Angles near 0 and pi have a singular gradient, so stay a little away
   * from the edge and check the value only where rounding matters. */
  for (t = 0; t < n; t++) {
    double th = M_PI-(drand48()+0.5)*1e-9;
    double loc[3][3] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                        {1.2*cos(th), 1.2*sin(th), 0.0}};
    double P[3][3], g[3][3], val;
    place(loc, P, 3);
    val = my_getangle(P[0], P[1], P[2], g[0], g[1], g[2]);
    if (!isfinite(val)) bad_angle++;
  }

  if (bad_dihed || bad_angle) {
    fprintf(stderr, "non-finite or wrong: %ld dihedrals, %ld angles\n",
            bad_dihed, bad_angle);
    return 1;
  }
  return 0;
}
