SHELL = /bin/sh

CC = gcc
AR = ar

LIB = libcparammodel.a
BIN_PARAMS = params
BIN_TEST = test

CCFLAGS = -c -g
CLFLAGS = -o
ARFLAGS = rcs

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

all: test.o sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o params.o
	$(AR) $(ARFLAGS) $(LIB) sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o
	$(CC) params.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_TEST)

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o $(BIN_PARAMS) $(BIN_TEST)
