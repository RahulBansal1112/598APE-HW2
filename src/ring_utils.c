#include "ring_utils.h"
#include "poly_utils.h"

Poly ring_add_mod(Poly* x, Poly* y, double modulus, Poly* poly_mod) {
  Poly sum = poly_add(x, y);

  Poly sum_mod = coeff_mod(&sum, modulus);
  // sum_mod = coeff_mod(sum_mod, modulus);

  Poly quot, rem;
  poly_divmod(&sum_mod, poly_mod, &quot, &rem);

  Poly rem_mod = coeff_mod(&rem, modulus);
  // rem_mod = coeff_mod(rem_mod, modulus);

  return rem_mod;
}

Poly ring_mul_mod(Poly* x, Poly* y, double modulus, Poly* poly_mod) {
  Poly prod = poly_mul(x, y);

  Poly prod_mod = coeff_mod(&prod, modulus);
  // prod_mod = coeff_mod(prod_mod, modulus);

  Poly quot, rem;
  poly_divmod(&prod_mod, poly_mod, &quot, &rem);

  Poly rem_mod = coeff_mod(&rem, modulus);
  // rem_mod = coeff_mod(rem_mod, modulus);

  return rem_mod;
}

Poly ring_mul_no_mod_q(Poly* x, Poly* y, Poly* poly_mod) {
  Poly prod = poly_mul(x, y);

  Poly quot, rem;
  poly_divmod(&prod, poly_mod, &quot, &rem);

  return rem;
}

Poly ring_add_no_mod_q(Poly* x, Poly* y, Poly* poly_mod) {
  Poly sum = poly_add(x, y);

  Poly quot, rem;
  poly_divmod(&sum, poly_mod, &quot, &rem);

  return rem;
}

Poly ring_mul_poly_mod(Poly* x, Poly* y, Poly* poly_mod) {
  Poly prod = poly_mul(x, y);
  Poly quot, rem;
  poly_divmod(&prod, poly_mod, &quot, &rem);
  return rem;
}

Poly ring_add_poly_mod(Poly* x, Poly* y, Poly* poly_mod) {
  Poly sum = poly_add(x, y);
  Poly quot, rem;
  poly_divmod(&sum, poly_mod, &quot, &rem);
  return rem;
}