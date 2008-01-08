SHELL = /bin/sh

CC = gcc

BIN_PARAMS = params
BIN_EXAMPLES = examples

CCFLAGS = -O2 -c -g
CLFLAGS = -o

%.o : %.c
	${CC} $(CCFLAGS) `pkg-config --cflags cparamlib` $< -o $@

all: examples.o params.o
	$(CC) params.c $(CLFLAGS) $(BIN_PARAMS) `pkg-config --libs --cflags cparamlib`
	$(CC) examples.c $(CLFLAGS) $(BIN_EXAMPLES) `pkg-config --libs --cflags cparamlib`

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.o $(BIN_PARAMS) $(BIN_EXAMPLES)
