#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <stdio.h>

Poly create_poly(void) {
  Poly p;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    p.coeffs[i] = 0.0;
  }

  return p;
}

double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(Poly* p) {
  for (int64_t i = MAX_POLY_DEGREE - 1; i >= 0; i--) {
    if (fabs(p->coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(Poly* p, int64_t degree) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p->coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  p->coeffs[degree] = value;
}

Poly coeff_mod(Poly* p, double modulus) {
  // Poly out = create_poly();
  Poly out;

  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(p->coeffs[i]) > 1e-9) {
      double rounded = round(p->coeffs[i]);
      double m = positive_fmod(rounded, modulus);
      out.coeffs[i] = m;
    } else {
      out.coeffs[i] = 0;
    }
  }
  return out;
}

Poly poly_add(Poly* a, Poly* b) {
  // Poly sum = create_poly();
  Poly sum;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    sum.coeffs[i] = a->coeffs[i] + b->coeffs[i];
  }
  return sum;
}

Poly poly_mul_scalar(Poly* p, double scalar) {
  // Poly res = create_poly();
  Poly res;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    res.coeffs[i] = p->coeffs[i] * scalar;
  }
  return res;
}

Poly poly_mul(Poly* a, Poly* b) {
  Poly res = create_poly();

  double deg_a = poly_degree(a);
  double deg_b = poly_degree(b);

  // #pragma omp parallel for
  for (int i = 0; i <= deg_a; i++) {
    if (fabs(a->coeffs[i]) > 1e-9) {
      double a_coeff = a->coeffs[i];
      for (int j = 0; j <= deg_b; j++) {
        if (fabs(b->coeffs[j]) > 1e-9) {
          assert(i + j < MAX_POLY_DEGREE);
          res.coeffs[i + j] += a_coeff * b->coeffs[j];
        }
      }
    }
  }
  return res;
}

void poly_divmod(Poly* num, Poly* den, Poly *quot, Poly *rem) {
  // In our case `den` should always be (x^n + 1)

  // assert(poly_degree(den) > 0 || fabs(get_coeff(den, 0)) > 1e-9);

  // size_t ndeg = poly_degree(num);
  // size_t ddeg = poly_degree(den);

  // *quot = create_poly();
  // *rem = *num;

  // if (ndeg < ddeg) {
  //   return;
  // }

  // double d_lead = get_coeff(den, ddeg);

  // assert(fabs(d_lead) > 1e-9);

  // for (int64_t k = ndeg - ddeg; k >= 0; --k) {
  //   int64_t target_deg = ddeg + k;
  //   // double r_coeff = get_coeff(*rem, target_deg);
  //   double r_coeff = rem->coeffs[target_deg];
  //   double coeff = trunc(round(r_coeff) / round(d_lead));
  //   quot->coeffs[k] += coeff;

  //   for (int i = 0; i <= ddeg; i++) {
  //     if (fabs(den->coeffs[i]) > 1e-9) {
  //       int64_t deg = i + k;
  //       assert(deg < MAX_POLY_DEGREE);
  //       rem->coeffs[deg] -= coeff * den->coeffs[i];
  //     }
  //   }
  // }

  // assert(poly_degree(rem) < poly_degree(den));

  *rem = *num;
  int n = poly_degree(den);
  for (int i = n; i < MAX_POLY_DEGREE; i++) {
    rem->coeffs[i - n] -= rem->coeffs[i];
    rem->coeffs[i] = 0;
  }
}

Poly poly_round_div_scalar(Poly* x, double divisor) {
  // Poly out = create_poly();
  Poly out;
  assert(fabs(divisor) > 1e-9);

  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double v = x->coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  return out;
}
