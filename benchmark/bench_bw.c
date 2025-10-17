#define STB_IMAGE_IMPLEMENTATION
#include "../external/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../external/stb_image_write.h"

#include "../src/he.h"
#include "../src/poly_utils.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define TILE_SIZE 1024

typedef struct {
  uint8_t *data;
  int width;
  int height;
  int channels;
} Image;

static Image load_image(const char *path) {
  Image img;
  img.data = stbi_load(path, &img.width, &img.height, &img.channels, 0);
  if (!img.data) {
    fprintf(stderr, "Failed to load image: %s\n", path);
    exit(1);
  }
  return img;
}

static void free_image(Image img) { stbi_image_free(img.data); }

static void save_image(const char *path, Image img) {
  stbi_write_png(path, img.width, img.height, img.channels, img.data,
                 img.width * img.channels);
}

static int64_t extended_gcd(int64_t a, int64_t b, int64_t *x, int64_t *y) {
  if (b == 0) {
    *x = 1;
    *y = 0;
    return a;
  }
  int64_t x1, y1;
  int64_t gcd = extended_gcd(b, a % b, &x1, &y1);
  *x = y1;
  *y = x1 - (a / b) * y1;
  return gcd;
}

static int64_t mod_inverse(int64_t a, int64_t m) {
  int64_t x, y;
  int64_t gcd = extended_gcd(a, m, &x, &y);
  if (gcd != 1) {
    return -1;
  }
  int64_t result = (x % m + m) % m;
  return result;
}

static void rgb_to_grayscale_plain(uint8_t *input, uint8_t *output, int width,
                                   int height, int channels) {
  for (int i = 0; i < width * height; i++) {
    if (channels >= 3) {
      uint8_t r = input[i * channels + 0];
      uint8_t g = input[i * channels + 1];
      uint8_t b = input[i * channels + 2];
      output[i] = (uint8_t)((r + g + b) / 3);
    } else {
      output[i] = input[i * channels];
    }
  }
}

static void rgb_to_grayscale_fhe(Ciphertext *r_enc, Ciphertext *g_enc,
                                 Ciphertext *b_enc, Ciphertext *output_enc,
                                 int total_pixels, int64_t q, int64_t t,
                                 Poly poly_mod) {
  int64_t inv3 = mod_inverse(3, t);
  assert(inv3 != -1 &&
         "3 has no modular inverse modulo t; choose t coprime with 3");
  for (int i = 0; i < total_pixels; i++) {
    Ciphertext sum = add_cipher(r_enc[i], g_enc[i], q, poly_mod);
    sum = add_cipher(sum, b_enc[i], q, poly_mod);
    output_enc[i] = mul_plain(sum, q, t, poly_mod, inv3);
  }
}

int main(int argc, char **argv) {
  srand(42);
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input_image>\n", argv[0]);
    return 1;
  }

  const char *input_path = argv[1];

  size_t n = 1u << 4;
  int64_t q = 1ll << 30;
  int64_t t = 769;

  printf("Loading image: %s\n", input_path);
  Image img = load_image(input_path);
  printf("Image size: %dx%d, channels: %d\n", img.width, img.height,
         img.channels);

  if (img.channels < 3) {
    fprintf(stderr, "Error: Image must have at least 3 channels (RGB)\n");
    free_image(img);
    return 1;
  }

  int total_pixels = img.width * img.height;

  Poly poly_mod = create_poly();
  set_coeff(&poly_mod, 0, 1);
  set_coeff(&poly_mod, n, 1);

  printf("Generating keys...\n");
  KeyPair keys = keygen(n, q, poly_mod);
  PublicKey pk = keys.pk;
  SecretKey sk = keys.sk;

  uint8_t *fhe_gray = (uint8_t *)malloc(total_pixels * sizeof(uint8_t));

  printf("Encrypting RGB channels...\n");
  int num_tiles = (total_pixels + TILE_SIZE - 1) / TILE_SIZE;

  double enc_start = omp_get_wtime();
  #pragma omp parallel for schedule(static)
  for (int tile = 0; tile < num_tiles; tile++) {
    int start = tile * TILE_SIZE;
    int end = start + TILE_SIZE;
    if (end > total_pixels) end = total_pixels;
    int tile_pixels = end - start;

    Ciphertext *r_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *g_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *b_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));

    for (int i = 0; i < tile_pixels; i++) {
      int idx = start + i;
      uint8_t r = img.data[idx * img.channels + 0];
      uint8_t g = img.data[idx * img.channels + 1];
      uint8_t b = img.data[idx * img.channels + 2];

      r_enc[i] = encrypt(pk, n, q, poly_mod, t, r);
      g_enc[i] = encrypt(pk, n, q, poly_mod, t, g);
      b_enc[i] = encrypt(pk, n, q, poly_mod, t, b);
    }

    free(r_enc);
    free(g_enc);
    free(b_enc);
  }
  double enc_end = omp_get_wtime();
  double enc_time = enc_end - enc_start;

  printf("Applying FHE grayscale conversion (R+G+B)/3...\n");
  double fhe_start = omp_get_wtime();
  #pragma omp parallel for schedule(static)
  for (int tile = 0; tile < num_tiles; tile++) {
    int start = tile * TILE_SIZE;
    int end = start + TILE_SIZE;
    if (end > total_pixels) end = total_pixels;
    int tile_pixels = end - start;

    Ciphertext *r_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *g_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *b_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *gray_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));

    for (int i = 0; i < tile_pixels; i++) {
      int idx = start + i;
      uint8_t r = img.data[idx * img.channels + 0];
      uint8_t g = img.data[idx * img.channels + 1];
      uint8_t b = img.data[idx * img.channels + 2];

      r_enc[i] = encrypt(pk, n, q, poly_mod, t, r);
      g_enc[i] = encrypt(pk, n, q, poly_mod, t, g);
      b_enc[i] = encrypt(pk, n, q, poly_mod, t, b);
    }

    rgb_to_grayscale_fhe(r_enc, g_enc, b_enc, gray_enc, tile_pixels, q, t,
                         poly_mod);

    free(r_enc);
    free(g_enc);
    free(b_enc);
    free(gray_enc);
  }
  double fhe_end = omp_get_wtime();
  double fhe_time = fhe_end - fhe_start;

  printf("Decrypting FHE grayscale result...\n");
  double dec_start = omp_get_wtime();
  #pragma omp parallel for schedule(static)
  for (int tile = 0; tile < num_tiles; tile++) {
    int start = tile * TILE_SIZE;
    int end = start + TILE_SIZE;
    if (end > total_pixels) end = total_pixels;
    int tile_pixels = end - start;

    Ciphertext *r_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *g_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *b_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));
    Ciphertext *gray_enc = (Ciphertext *)malloc(tile_pixels * sizeof(Ciphertext));

    for (int i = 0; i < tile_pixels; i++) {
      int idx = start + i;
      uint8_t r = img.data[idx * img.channels + 0];
      uint8_t g = img.data[idx * img.channels + 1];
      uint8_t b = img.data[idx * img.channels + 2];

      r_enc[i] = encrypt(pk, n, q, poly_mod, t, r);
      g_enc[i] = encrypt(pk, n, q, poly_mod, t, g);
      b_enc[i] = encrypt(pk, n, q, poly_mod, t, b);
    }

    rgb_to_grayscale_fhe(r_enc, g_enc, b_enc, gray_enc, tile_pixels, q, t,
                         poly_mod);

    int64_t th1 = (t + 2) / 3;
    int64_t th2 = (2 * t + 2) / 3;
    for (int i = 0; i < tile_pixels; i++) {
      int idx = start + i;
      int64_t val = decrypt(sk, n, q, poly_mod, t, gray_enc[i]);
      if (val >= th2) {
        val -= th2;
      } else if (val >= th1) {
        val -= th1;
      }
      if (val > 255)
        val = 255;
      if (val < 0)
        val = 0;
      fhe_gray[idx] = (uint8_t)val;
    }

    free(r_enc);
    free(g_enc);
    free(b_enc);
    free(gray_enc);
  }
  double dec_end = omp_get_wtime();
  double dec_time = dec_end - dec_start;

  printf("Computing plaintext reference grayscale...\n");
  uint8_t *plain_gray = (uint8_t *)malloc(total_pixels * sizeof(uint8_t));
  clock_t plain_start = clock();
  rgb_to_grayscale_plain(img.data, plain_gray, img.width, img.height,
                         img.channels);
  clock_t plain_end = clock();
  double plain_time = ((double)(plain_end - plain_start)) / CLOCKS_PER_SEC;

  printf("Computing L2 norm error...\n");
  double l2_error = 0.0;
  int max_diff = 0;
  int num_errors = 0;
  for (int i = 0; i < total_pixels; i++) {
    int diff = abs((int)fhe_gray[i] - (int)plain_gray[i]);
    if (diff > 0)
      num_errors++;
    if (diff > max_diff)
      max_diff = diff;
    l2_error += (double)diff * (double)diff;
  }
  l2_error = sqrt(l2_error);
  double avg_error = l2_error / total_pixels;

  printf("First 10 pixels comparison:\n");
  for (int i = 0; i < 10 && i < total_pixels; i++) {
    printf("  Pixel %d: plain=%d, fhe=%d, diff=%d\n", i, plain_gray[i],
           fhe_gray[i], abs((int)fhe_gray[i] - (int)plain_gray[i]));
  }

  printf("\n=== Results ===\n");
  printf("Encryption time: %.4f s (%.2f ms/pixel)\n", enc_time,
         enc_time * 1000.0 / total_pixels);
  printf("FHE grayscale conversion time: %.4f s\n", fhe_time);
  printf("Decryption time: %.4f s (%.2f ms/pixel)\n", dec_time,
         dec_time * 1000.0 / total_pixels);
  printf("Plaintext grayscale time: %.4f s\n", plain_time);
  printf("L2 norm error: %.4f\n", l2_error);
  printf("Average error per pixel: %.4f\n", avg_error);
  printf("Max pixel difference: %d\n", max_diff);
  printf("Pixels with errors: %d/%d (%.1f%%)\n", num_errors, total_pixels,
         100.0 * num_errors / total_pixels);

  save_image("output/original.png",
             (Image){img.data, img.width, img.height, img.channels});
  save_image("output/fhe_grayscale.png",
             (Image){fhe_gray, img.width, img.height, 1});
  save_image("output/plaintext_grayscale.png",
             (Image){plain_gray, img.width, img.height, 1});

  printf("\nSaved outputs:\n");
  printf("  output/original.png             (original RGB input)\n");
  printf("  output/fhe_grayscale.png        (FHE encrypted RGB -> grayscale -> "
         "decrypted)\n");
  printf("  output/plaintext_grayscale.png  (plaintext reference grayscale)\n");

  free(fhe_gray);
  free(plain_gray);
  free_image(img);

  return 0;
}