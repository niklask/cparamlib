SHELL = /bin/sh

CC = gcc
AR = ar

LIB = libcparammodel.a

CCFLAGS = -c -g
CLFLAGS = -o
ARFLAGS = rcs

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

all: test.o flux.o gamma.o params.o
	$(AR) $(ARFLAGS) $(LIB) flux.o gamma.o
	$(CC) params.c -L. -lm -lcparammodel $(CLFLAGS) params
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) test

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o
