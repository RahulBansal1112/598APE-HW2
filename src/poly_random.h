#ifndef POLY_RANDOM_H
#define POLY_RANDOM_H

#include "types.h"
#include <stdint.h>

Poly gen_binary_poly(size_t size);

Poly gen_uniform_poly(size_t size, double modulus);

Poly gen_normal_poly(size_t size, double mean, double stddev);

#endif