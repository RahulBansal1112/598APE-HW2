#ifndef HE_H
#define HE_H

#include "types.h"
#include <stdint.h>

typedef struct {
  PublicKey pk;
  SecretKey sk;
} KeyPair;

KeyPair keygen(size_t n, double q, Poly poly_mod);

Ciphertext encrypt(PublicKey pk, size_t n, double q, Poly poly_mod, double t,
                   double pt);

double decrypt(SecretKey sk, size_t n, double q, Poly poly_mod, double t,
               Ciphertext ct);

Poly encode_plain_integer(double t, double pt);

Ciphertext add_plain(Ciphertext ct, double q, double t, Poly poly_mod,
                     double pt);

Ciphertext add_cipher(Ciphertext c1, Ciphertext c2, double q, Poly poly_mod);

Ciphertext mul_plain(Ciphertext ct, double q, double t, Poly poly_mod,
                     double pt);

EvalKey evaluate_keygen(SecretKey sk, size_t n, double q, Poly poly_mod,
                        double p);

Ciphertext mul_cipher(Ciphertext c1, Ciphertext c2, double q, double t,
                      double p, Poly poly_mod, EvalKey rlk);

#endif