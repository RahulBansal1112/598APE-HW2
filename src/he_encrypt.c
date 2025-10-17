#include "he.h"
#include "poly_random.h"
#include "poly_utils.h"
#include "ring_utils.h"

Poly encode_plain_integer(double t, double pt) {
  Poly m = create_poly();
  double v = positive_fmod(pt, t);
  m.coeffs[0] = v;
  return m;
}

Ciphertext encrypt(PublicKey pk, size_t n, double q, Poly poly_mod, double t,
                   double pt) {
  Poly u = gen_binary_poly(n);
  Poly e1 = gen_normal_poly(n, 0.0, 1.0);
  Poly e2 = gen_normal_poly(n, 0.0, 1.0);

  Poly bu = ring_mul_mod(&pk.b, &u, q, &poly_mod);
  Poly au = ring_mul_mod(&pk.a, &u, q, &poly_mod);

  Poly c0 = ring_add_mod(&bu, &e1, q, &poly_mod);
  Poly c1 = ring_add_mod(&au, &e2, q, &poly_mod);

  double scaled_pt = floor(q / t) * positive_fmod(pt, t);
  c0.coeffs[0] = positive_fmod(c0.coeffs[0] + scaled_pt, q);

  Ciphertext ct;
  ct.c0 = c0;
  ct.c1 = c1;
  return ct;
}