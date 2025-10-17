#ifndef RING_UTILS_H
#define RING_UTILS_H

#include "types.h"
#include <stdint.h>

Poly ring_add_mod(Poly* x, Poly* y, double modulus, Poly* poly_mod);

Poly ring_mul_mod(Poly* x, Poly* y, double modulus, Poly* poly_mod);

Poly ring_mul_no_mod_q(Poly* x, Poly* y, Poly* poly_mod);

Poly ring_add_no_mod_q(Poly* x, Poly* y, Poly* poly_mod);

Poly ring_add_poly_mod(Poly* x, Poly* y, Poly* poly_mod);

Poly ring_mul_poly_mod(Poly* x, Poly* y, Poly* poly_mod);

#endif