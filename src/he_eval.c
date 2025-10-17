#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <assert.h>
#include <math.h>

Ciphertext add_plain(Ciphertext ct, double q, double t, Poly poly_mod,
                     double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly scaled_m = poly_mul_scalar(&m, q / t);
  Poly new_c0 = ring_add_mod(&ct.c0, &scaled_m, q, &poly_mod);

  Ciphertext result;
  result.c0 = new_c0;
  result.c1 = ct.c1;
  return result;
}

Ciphertext add_cipher(Ciphertext c1, Ciphertext c2, double q, Poly poly_mod) {
  Poly new_c0 = ring_add_mod(&c1.c0, &c2.c0, q, &poly_mod);
  Poly new_c1 = ring_add_mod(&c1.c1, &c2.c1, q, &poly_mod);

  Ciphertext result;
  result.c0 = new_c0;
  result.c1 = new_c1;
  return result;
}

Ciphertext mul_plain(Ciphertext ct, double q, double t, Poly poly_mod,
                     double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly new_c0 = ring_mul_mod(&ct.c0, &m, q, &poly_mod);
  Poly new_c1 = ring_mul_mod(&ct.c1, &m, q, &poly_mod);

  Ciphertext result;
  result.c0 = new_c0;
  result.c1 = new_c1;
  return result;
}

Ciphertext mul_cipher(Ciphertext c1, Ciphertext c2, double q, double t,
                      double p, Poly poly_mod, EvalKey rlk) {
  Poly c0_prod = ring_mul_no_mod_q(&c1.c0, &c2.c0, &poly_mod);
  Poly c1_left = ring_mul_no_mod_q(&c1.c0, &c2.c1, &poly_mod);
  Poly c1_right = ring_mul_no_mod_q(&c1.c1, &c2.c0, &poly_mod);
  Poly c1_sum = ring_add_no_mod_q(&c1_left, &c1_right, &poly_mod);
  Poly c2_prod = ring_mul_no_mod_q(&c1.c1, &c2.c1, &poly_mod);

  Poly c0_res = create_poly();
  Poly c1_res = create_poly();
  Poly c2_res = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(c0_prod.coeffs[i]) > 1e-9) {
      c0_res.coeffs[i] = round(t * c0_prod.coeffs[i] / q);
    }
    if (fabs(c1_sum.coeffs[i]) > 1e-9) {
      c1_res.coeffs[i] = round(t * c1_sum.coeffs[i] / q);
    }
    if (fabs(c2_prod.coeffs[i]) > 1e-9) {
      c2_res.coeffs[i] = round(t * c2_prod.coeffs[i] / q);
    }
  }

  Poly c0_modq = coeff_mod(&c0_res, q);
  Poly c1_modq = coeff_mod(&c1_res, q);
  Poly c2_modq = coeff_mod(&c2_res, q);

  // Relinearization
  Poly prod_b = ring_mul_no_mod_q(&rlk.b, &c2_modq, &poly_mod);
  Poly prod_a = ring_mul_no_mod_q(&rlk.a, &c2_modq, &poly_mod);

  Poly div_b = create_poly();
  Poly div_a = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double vb = prod_b.coeffs[i];
    double va = prod_a.coeffs[i];
    if (fabs(vb) > 1e-9) {
      div_b.coeffs[i] = round(vb / p);
    }
    if (fabs(va) > 1e-9) {
      div_a.coeffs[i] = round(va / p);
    }
  }

  Poly c20_modq = coeff_mod(&div_b, q);
  Poly c21_modq = coeff_mod(&div_a, q);

  Poly new_c0 = ring_add_mod(&c0_modq, &c20_modq, q, &poly_mod);
  Poly new_c1 = ring_add_mod(&c1_modq, &c21_modq, q, &poly_mod);

  Ciphertext out;
  out.c0 = new_c0;
  out.c1 = new_c1;
  return out;
}