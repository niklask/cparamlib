SHELL = /bin/sh

CC = gcc

BIN_PARAMS = parameters
BIN_EX1 = example1
BIN_SPECTRUM = spectrum

CCFLAGS = -O2 -c -g
CLFLAGS = -o

%.o : %.c
	${CC} $(CCFLAGS) `pkg-config --cflags cparamlib` $< -o $@

all: example1.o parameters.o spectrum.o
	$(CC) parameters.c $(CLFLAGS) $(BIN_PARAMS) `pkg-config --libs --cflags cparamlib`
	$(CC) example1.c $(CLFLAGS) $(BIN_EX1) `pkg-config --libs --cflags cparamlib`
	$(CC) spectrum.c $(CLFLAGS) $(BIN_SPECTRUM) `pkg-config --libs --cflags cparamlib`

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.o *.csv $(BIN_PARAMS) $(BIN_EXAMPLES) $(BIN_SPECTRUM)
