#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "types.h"
#include <math.h>
#include <stdint.h>
#include <assert.h>

double positive_fmod(double x, double m);
// __attribute__((always_inline)) inline double positive_fmod(double x, double m) {
//   assert(m > 0.0);
//   double r = fmod(x, m);
//   if (r < 0.0)
//     r += m;
//   return r;
// }

int64_t poly_degree(Poly *p);

double get_coeff(Poly *p, int64_t degree);

void set_coeff(Poly *p, int64_t degree, double value);

Poly coeff_mod(Poly* p, double modulus);

Poly poly_add(Poly* a, Poly* b);

Poly poly_mul_scalar(Poly* p, double scalar);

Poly poly_mul(Poly* a, Poly* b);

void poly_divmod(Poly* numerator, Poly* denominator, Poly *quotient,
                 Poly *remainder);

Poly create_poly(void);

#endif