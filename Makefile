CC := gcc
CFLAGS := -O3 -lm -g -Werror -std=c99 -fopenmp -DNDEBUG
OBJ_DIR := ./bin/

C_FILES := $(wildcard src/*.c)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(C_FILES:.c=.o)))

all: $(OBJ_DIR) $(OBJ_FILES)
	$(CC) ./main.c -o ./main.exe $(OBJ_FILES) $(CFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o: src/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

bench_matmul.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_matmul.c
	$(CC) ./benchmark/bench_matmul.c -o $@ $(OBJ_FILES) $(CFLAGS)

bench_bw.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_bw.c
	$(CC) ./benchmark/bench_bw.c -o $@ $(OBJ_FILES) $(CFLAGS)

bench_sobel.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_sobel.c
	$(CC) ./benchmark/bench_sobel.c -o $@ $(OBJ_FILES) $(CFLAGS)

clean:
	rm -f ./*.exe
	rm -f ./*.o
	rm -rf $(OBJ_DIR)
