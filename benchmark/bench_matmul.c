#include "../src/he.h"
#include "../src/poly_utils.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

static int64_t **alloc_matrix(size_t rows, size_t cols) {
  int64_t **M = (int64_t **)malloc(rows * sizeof(int64_t *));
  for (size_t i = 0; i < rows; i++) {
    M[i] = (int64_t *)calloc(cols, sizeof(int64_t));
  }
  return M;
}

static void free_matrix(int64_t **M, size_t rows) {
  for (size_t i = 0; i < rows; i++) {
    free(M[i]);
  }
  free(M);
}

static Ciphertext **alloc_ct_matrix(size_t rows, size_t cols) {
  Ciphertext **M = (Ciphertext **)malloc(rows * sizeof(Ciphertext *));
  for (size_t i = 0; i < rows; i++) {
    M[i] = (Ciphertext *)malloc(cols * sizeof(Ciphertext));
  }
  return M;
}

static void free_ct_matrix(Ciphertext **M, size_t rows) {
  for (size_t i = 0; i < rows; i++) {
    free(M[i]);
  }
  free(M);
}

int main(int argc, char **argv) {
  srand(42);
  int mode = 1; // 0 is ct * pt mode, 1 is ct * ct mode

  // Please report runtimes on the following parameters
  size_t dim = 32;
  size_t n = 1u << 4;
  int64_t q = 1ll << 32; // Large q like 1<<40 can make results unreliable
  int64_t t = 1ll << 8;

  if (argc >= 2)
    mode = atoi(argv[1]);
  if (argc >= 3)
    dim = (size_t)strtoull(argv[2], NULL, 10);
  if (argc >= 4)
    n = (size_t)strtoull(argv[3], NULL, 10);

  printf("Matrix size: %zux%zu, Mode: %d (%s)\n", dim, dim, mode,
         mode == 0 ? "ct*pt" : "ct*ct");

  Poly poly_mod = create_poly();
  set_coeff(&poly_mod, 0, 1.0);
  set_coeff(&poly_mod, n, 1.0);

  KeyPair keys = keygen(n, q, poly_mod);
  PublicKey pk = keys.pk;
  SecretKey sk = keys.sk;

  // Generate plaintext matrices A, B
  int64_t **A = alloc_matrix(dim, dim);
  int64_t **B = alloc_matrix(dim, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      A[i][j] = rand() % t;
      B[i][j] = rand() % t;
    }
  }

  // Plaintext reference C = A * B (mod t)
  int64_t **C_ref = alloc_matrix(dim, dim);
  double ref_start = omp_get_wtime();

  // Plaintext matrix multiplication
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < dim; ++i) {
    for (size_t k = 0; k < dim; ++k) {
      int64_t acc = 0;
      for (size_t j = 0; j < dim; ++j) {
        acc += A[i][j] * B[j][k];
      }
      C_ref[i][k] = acc % t;
    }
  }
  double ref_end = omp_get_wtime();
  double ref_sec = ref_end - ref_start;

  // Start timing the entire encrypted workflow
  double enc_start = omp_get_wtime();

  // Encrypt B (and optionally A)
  Ciphertext **B_enc = alloc_ct_matrix(dim, dim);
  #pragma omp parallel for schedule(static)
  for (size_t j = 0; j < dim; ++j) {
    for (size_t k = 0; k < dim; ++k) {
      B_enc[j][k] = encrypt(pk, n, q, poly_mod, t, B[j][k]);
    }
  }

  Ciphertext **A_enc = NULL;
  EvalKey evk;
  double p = pow(q, 2.0);
  if (mode == 1) {
    A_enc = alloc_ct_matrix(dim, dim);
    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        A_enc[i][j] = encrypt(pk, n, q, poly_mod, t, A[i][j]);
      }
    }
    evk = evaluate_keygen(sk, n, q, poly_mod, p);
  }

  // Encrypted matmul
  Ciphertext **C_enc = alloc_ct_matrix(dim, dim);

  if (mode == 0) {
    // Mode 0: ct * pt matmul: C_enc[i][k] = sum_j A[i][j] * Enc(B[j][k])
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dim; ++i) {
      for (size_t k = 0; k < dim; ++k) {
        int first = 1;
        Ciphertext acc_ct;

        for (size_t j = 0; j < dim; ++j) {
          Ciphertext term = mul_plain(B_enc[j][k], q, t, poly_mod, A[i][j]);
          if (first) {
            acc_ct = term;
            first = 0;
          } else {
            acc_ct = add_cipher(acc_ct, term, q, poly_mod);
          }
        }
        C_enc[i][k] = acc_ct;
      }
    }
  } else {
    // Mode 1: ct * ct matmul: C_enc[i][k] = sum_j Enc(A[i][j]) * Enc(B[j][k])
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dim; ++i) {
      for (size_t k = 0; k < dim; ++k) {
        int first = 1;
        Ciphertext acc_ct;

        for (size_t j = 0; j < dim; ++j) {
          Ciphertext term =
              mul_cipher(A_enc[i][j], B_enc[j][k], q, t, p, poly_mod, evk);
          if (first) {
            acc_ct = term;
            first = 0;
          } else {
            acc_ct = add_cipher(acc_ct, term, q, poly_mod);
          }
        }
        C_enc[i][k] = acc_ct;
      }
    }
  }

  double enc_end = omp_get_wtime();
  double enc_sec = enc_end - enc_start;

  // Decrypt result matrix
  int64_t **C_dec = alloc_matrix(dim, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t k = 0; k < dim; ++k) {
      C_dec[i][k] = decrypt(sk, n, q, poly_mod, t, C_enc[i][k]);
    }
  }

  // Relative error (Frobenius): ||C_dec - C_ref||_F / ||C_ref||_F
  long double diff_acc = 0.0L;
  long double ref_acc = 0.0L;

  for (size_t i = 0; i < dim; ++i) {
    for (size_t k = 0; k < dim; ++k) {
      long double d = (long double)(C_dec[i][k]) - (long double)(C_ref[i][k]);
      diff_acc += d * d;
      long double rv = (long double)(C_ref[i][k]);
      ref_acc += rv * rv;
    }
  }

  double rel_err = 0.0;
  if (ref_acc > 0.0L) {
    rel_err = sqrt((double)(diff_acc / ref_acc));
  } else {
    rel_err = (diff_acc == 0.0L) ? 0.0 : DBL_MAX;
  }

  printf("ref_time_sec=%f, enc_time_sec=%f, rel_err=%f\n", ref_sec, enc_sec,
         rel_err);

  size_t show = (3 < dim) ? 3 : dim;
  printf("C_ref (top %zux%zu):\n", show, show);
  for (size_t i = 0; i < show; ++i) {
    for (size_t j = 0; j < show; ++j) {
      if (j)
        printf(" ");
      printf("%ld", C_ref[i][j]);
    }
    printf("\n");
  }

  printf("C_dec (top %zux%zu):\n", show, show);
  for (size_t i = 0; i < show; ++i) {
    for (size_t j = 0; j < show; ++j) {
      if (j)
        printf(" ");
      printf("%ld", C_dec[i][j]);
    }
    printf("\n");
  }

  free_matrix(A, dim);
  free_matrix(B, dim);
  free_matrix(C_ref, dim);
  free_matrix(C_dec, dim);
  free_ct_matrix(B_enc, dim);
  free_ct_matrix(C_enc, dim);
  if (mode == 1) {
    free_ct_matrix(A_enc, dim);
  }

  return 0;
}