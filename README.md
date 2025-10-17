# 598APE-HW2: Homomorphic Encryption

This is a C implementation of a SHE scheme based on the [Fan-Vercauteren](https://eprint.iacr.org/2012/144) RLWE-based approach and [toy implementations in Python](https://bit-ml.github.io/blog/post/homomorphic-encryption-toy-implementation-in-python/). Your task is to build and run it, then consider performance improvements.

Warning: This code is for instruction only. Do not use this implementation in any practical setting!

To compile the program run:
```bash
make -j
```

To clean existing build artifacts run:
```bash
make clean
```

Once built, run:
```bash
./main.exe  # A short demonstration
```

## Benchmarks

To build and run the benchmarks:

```bash
make bench_matmul.exe
./bench_matmul.exe 0 16  # 16x16 Ciphertext * Plaintext
./bench_matmul.exe 1 32  # 32x32 Ciphertext * Ciphertext

# Build and run B+W image converter benchmark
make bench_bw.exe
./bench_bw.exe inputs/bird.jpg

# Build and run Sobel filter benchmark
make bench_sobel.exe
./bench_sobel.exe inputs/bird.jpg
```

## Docker

For ease of use and installation, we provide a docker image capable of running and building code here. The source docker file is in /docker (which is essentially a list of commands to build an OS state from scratch). It contains the dependent compilers, and some other nice things.

You can build this yourself manually by running `cd docker && docker build -t <myusername>/598ape`. Alternatively we have pushed a pre-built version to `wsmoses/598ape` on Dockerhub.

You can then use the Docker container to build and run your code. If you run `./dockerrun.sh` you will enter an interactive bash session with all the packages from docker installed (that script by default uses `wsmoses/598ape`, feel free to replace it with whatever location you like if you built from scratch). The current directory (aka this folder) is mounted within `/host`. Any files you create on your personal machine will be available there, and anything you make in the container in that folder will be available on your personal machine.

## References

- Fan, J., & Vercauteren, F. (2012). *Somewhat Practical Fully Homomorphic Encryption*. Cryptology ePrint Archive, Paper 2012/144. https://eprint.iacr.org/2012/144

- *Homomorphic Encryption: Toy Implementation in Python*. https://bit-ml.github.io/blog/post/homomorphic-encryption-toy-implementation-in-python/

- Benaissa, A. *Build an HE Scheme from Scratch in Python*. https://www.ayoub-benaissa.com/blog/build-he-scheme-from-scratch-python
