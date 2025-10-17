#include "src/he.h"
#include "src/poly_utils.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
  srand(time(NULL));
  size_t n = 1u << 4;
  int64_t q = 1ll << 28;
  int64_t t = 1ll << 8;

  // poly_mod = X^n + 1
  Poly poly_mod = create_poly();
  set_coeff(&poly_mod, 0, 1.0);
  set_coeff(&poly_mod, n, 1.0);

  KeyPair keys = keygen(n, q, poly_mod);
  PublicKey pk = keys.pk;
  SecretKey sk = keys.sk;

  printf("[+] Public Key:\n\n");
  printf("\t pk.b: [");
  int first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(pk.b.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, pk.b.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n");

  printf("\t pk.a: [");
  first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(pk.a.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, pk.a.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n\n");

  int64_t pt1 = 73;
  int64_t pt2 = 20;
  int64_t cst1 = 7;
  int64_t cst2 = 5;

  Ciphertext ct1 = encrypt(pk, n, q, poly_mod, t, pt1);
  Ciphertext ct2 = encrypt(pk, n, q, poly_mod, t, pt2);

  printf("[+] Ciphertext ct1(%ld):\n\n", pt1);
  printf("\t ct1_0: [");
  first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(ct1.c0.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, ct1.c0.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n");

  printf("\t ct1_1: [");
  first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(ct1.c1.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, ct1.c1.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n\n");

  printf("[+] Ciphertext ct2(%ld):\n\n", pt2);
  printf("\t ct2_0: [");
  first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(ct2.c0.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, ct2.c0.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n");

  printf("\t ct2_1: [");
  first = 1;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(ct2.c1.coeffs[i]) > 1e-9) {
      if (!first)
        printf(", ");
      printf("%d:%.0f", i, ct2.c1.coeffs[i]);
      first = 0;
    }
  }
  printf("]\n\n");

  Ciphertext ct3 = add_plain(ct1, q, t, poly_mod, cst1);
  Ciphertext ct4 = mul_plain(ct2, q, t, poly_mod, cst2);
  Ciphertext ct5 = add_cipher(ct3, ct4, q, poly_mod);

  int64_t d3 = decrypt(sk, n, q, poly_mod, t, ct3);
  int64_t d4 = decrypt(sk, n, q, poly_mod, t, ct4);
  int64_t d5 = decrypt(sk, n, q, poly_mod, t, ct5);

  printf("[+] Decrypted ct3(ct1 + %ld): %ld\n", cst1, d3);
  printf("[+] Decrypted ct4(ct2 * %ld): %ld\n", cst2, d4);
  printf("[+] Decrypted ct5(ct1 + %ld + %ld * ct2): %ld\n", cst1, cst2, d5);

  int64_t expected = ((pt1 % t) * (pt2 % t)) % t;
  int64_t p = q * q;
  EvalKey rlk = evaluate_keygen(sk, n, q, poly_mod, p);
  Ciphertext ct7 = mul_cipher(ct1, ct2, q, t, p, poly_mod, rlk);
  int64_t d7 = decrypt(sk, n, q, poly_mod, t, ct7);
  printf("[+] Decrypted ct7(relin_v2 ct1*ct2): %ld (expected %ld)\n", d7,
         expected);
  if (d7 == expected) {
    printf("[OK] relin_v2 ct1 * ct2 mod t matches expected.\n");
  } else {
    printf("[FAIL] relin_v2 ct1 * ct2 mismatch.\n");
  }

  return 0;
}