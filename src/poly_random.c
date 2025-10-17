#include "poly_random.h"
#include "poly_utils.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

Poly gen_binary_poly(size_t size) {
  Poly p = create_poly();

  for (size_t i = 0; i < size && i < MAX_POLY_DEGREE; ++i) {
    double v = (rand() % 2) ? 1.0 : 0.0;
    p.coeffs[i] = v;
  }
  return p;
}

// Box-Muller transform
static double gen_normal(double mean, double stddev) {
  double u, v, s;
  do {
    u = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
    v = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
    s = u * u + v * v;
  } while (s >= 1.0 || s == 0.0);

  s = sqrt(-2.0 * log(s) / s);
  return mean + u * s * stddev;
}

Poly gen_normal_poly(size_t size, double mean, double stddev) {
  Poly p = create_poly();

  for (size_t i = 0; i < size && i < MAX_POLY_DEGREE; ++i) {
    double v = round(gen_normal(mean, stddev));
    p.coeffs[i] = v;
  }
  return p;
}

Poly gen_uniform_poly(size_t size, double modulus) {
  Poly p = create_poly();

  for (size_t i = 0; i < size && i < MAX_POLY_DEGREE; ++i) {
    double v = ((double)rand() / RAND_MAX) * modulus;
    p.coeffs[i] = v;
  }
  return p;
}