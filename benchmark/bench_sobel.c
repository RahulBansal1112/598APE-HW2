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

#define TILE_SIZE 64

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

static const int sobel_gx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
static const int sobel_gy[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

static void rgb_to_grayscale(uint8_t *input, uint8_t *output, int width,
                             int height, int channels) {
  #pragma omp parallel for schedule(static)
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

static void sobel_plain(uint8_t *input, uint8_t *output, int width,
                        int height) {
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      int gx = 0, gy = 0;

      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          int pixel = input[(y + ky) * width + (x + kx)];
          gx += pixel * sobel_gx[ky + 1][kx + 1];
          gy += pixel * sobel_gy[ky + 1][kx + 1];
        }
      }

      int magnitude = sqrt(gx * gx + gy * gy);
      if (magnitude > 255)
        magnitude = 255;
      output[y * width + x] = (uint8_t)magnitude;
    }
  }

  for (int x = 0; x < width; x++) {
    output[x] = 0;
    output[(height - 1) * width + x] = 0;
  }
  for (int y = 0; y < height; y++) {
    output[y * width] = 0;
    output[y * width + (width - 1)] = 0;
  }
}

Ciphertext encode_zero(int64_t q, Poly poly_mod) {
  Ciphertext ct;
  ct.c0 = encode_plain_integer(q, 0);
  ct.c1 = encode_plain_integer(q, 0);
  return ct;
}

static void sobel_fhe_tiled(uint8_t *gray, uint8_t *output, int width,
                              int height, PublicKey pk, SecretKey sk,
                              size_t n, int64_t q, int64_t t, Poly poly_mod,
                              double *enc_time, double *fhe_time, double *dec_time) {
  
  *enc_time = 0.0;
  *fhe_time = 0.0;
  *dec_time = 0.0;
  
  for (int x = 0; x < width; x++) {
    output[x] = 0;
    output[(height - 1) * width + x] = 0;
  }
  for (int y = 0; y < height; y++) {
    output[y * width] = 0;
    output[y * width + (width - 1)] = 0;
  }
  
  // #pragma omp parallel for schedule(static)
  for (int chunk_start = 1; chunk_start < height - 1; chunk_start += TILE_SIZE) {
    int chunk_end = chunk_start + TILE_SIZE;
    if (chunk_end > height - 1) chunk_end = height - 1;
    
    int row_start = chunk_start - 1;
    int row_end = chunk_end + 1;
    if (row_end > height) row_end = height;
    
    int chunk_height = row_end - row_start;
    int chunk_pixels = width * chunk_height;
    
    Ciphertext *chunk_enc = (Ciphertext *)malloc(chunk_pixels * sizeof(Ciphertext));
    
    // Encrypt chunk
    double enc_start_time = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < chunk_pixels; i++) {
      int global_idx = (row_start * width) + i;
      chunk_enc[i] = encrypt(pk, n, q, poly_mod, t, gray[global_idx]);
    }
    *enc_time += omp_get_wtime() - enc_start_time;
    
    // Apply Sobel to this chunk
    double fhe_start_time = omp_get_wtime();
    
    #pragma omp parallel for schedule(static)
    for (int y = chunk_start; y < chunk_end; y++) {
      for (int x = 1; x < width - 1; x++) {
        Ciphertext gx = encode_zero(q, poly_mod);
        Ciphertext gy = encode_zero(q, poly_mod);

        for (int ky = -1; ky <= 1; ky++) {
          for (int kx = -1; kx <= 1; kx++) {
            int global_y = y + ky;
            int local_y = global_y - row_start;
            int local_idx = local_y * width + (x + kx);
            
            Ciphertext pixel = chunk_enc[local_idx];

            int coeff_gx = sobel_gx[ky + 1][kx + 1];
            int coeff_gy = sobel_gy[ky + 1][kx + 1];

            if (coeff_gx != 0) {
              Ciphertext term = mul_plain(pixel, q, t, poly_mod, coeff_gx);
              gx = add_cipher(gx, term, q, poly_mod);
            }

            if (coeff_gy != 0) {
              Ciphertext term = mul_plain(pixel, q, t, poly_mod, coeff_gy);
              gy = add_cipher(gy, term, q, poly_mod);
            }
          }
        }

        // Decrypt immediately and store result
        Ciphertext result = add_cipher(gx, gy, q, poly_mod);
        int64_t val = decrypt(sk, n, q, poly_mod, t, result);
        if (val > t / 2)
          val = t - val;
        if (val > 255)
          val = 255;
        output[y * width + x] = (uint8_t)val;
      }
    }
    
    *fhe_time += omp_get_wtime() - fhe_start_time;
    
    free(chunk_enc);
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
  int64_t t = 1ll << 10;

  printf("Loading image: %s\n", input_path);
  Image img = load_image(input_path);
  printf("Image size: %dx%d, channels: %d\n", img.width, img.height,
         img.channels);

  uint8_t *gray = (uint8_t *)malloc(img.width * img.height * sizeof(uint8_t));
  rgb_to_grayscale(img.data, gray, img.width, img.height, img.channels);

  int total_pixels = img.width * img.height;

  Poly poly_mod = create_poly();
  set_coeff(&poly_mod, 0, 1);
  set_coeff(&poly_mod, n, 1);

  printf("Generating keys...\n");
  KeyPair keys = keygen(n, q, poly_mod);
  PublicKey pk = keys.pk;
  SecretKey sk = keys.sk;

  printf("Applying FHE Sobel edge detection (tiled)...\n");
  uint8_t *fhe_sobel = (uint8_t *)malloc(total_pixels * sizeof(uint8_t));
  
  double enc_time, fhe_time, dec_time;
  sobel_fhe_tiled(gray, fhe_sobel, img.width, img.height, pk, sk,
                    n, q, t, poly_mod, &enc_time, &fhe_time, &dec_time);

  printf("Computing plaintext Sobel edge detection...\n");
  uint8_t *plain_sobel = (uint8_t *)calloc(total_pixels, sizeof(uint8_t));
  double plain_start = omp_get_wtime();
  sobel_plain(gray, plain_sobel, img.width, img.height);
  double plain_end = omp_get_wtime();
  double plain_time = plain_end - plain_start;

  printf("\n=== Results ===\n");
  printf("Encryption time: %.4f s (%.2f ms/pixel)\n", enc_time,
         enc_time * 1000.0 / total_pixels);
  printf("FHE Sobel time: %.4f s\n", fhe_time);
  printf("Decryption time: %.4f s (%.2f ms/pixel)\n", dec_time,
         dec_time * 1000.0 / total_pixels);
  printf("Plaintext Sobel time: %.4f s\n", plain_time);

  save_image("output/sobel_fhe.png",
             (Image){fhe_sobel, img.width, img.height, 1});
  save_image("output/sobel_plain.png",
             (Image){plain_sobel, img.width, img.height, 1});
  save_image("output/grayscale.png", (Image){gray, img.width, img.height, 1});

  printf("\nSaved outputs:\n");
  printf("  output/grayscale.png     (grayscale input)\n");
  printf("  output/sobel_fhe.png     (FHE Sobel edges)\n");
  printf("  output/sobel_plain.png   (plaintext Sobel edges)\n");

  free(fhe_sobel);
  free(plain_sobel);
  free(gray);
  free_image(img);

  return 0;
}