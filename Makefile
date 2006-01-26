SHELL = /bin/sh

CC = gcc
AR = ar

LIB = libcparammodel.a

CCFLAGS = -c -g
CLFLAGS = -o
ARFLAGS = rcs

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

all: test.o flux.o gamma.o
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) test
	$(AR) $(ARFLAGS) $(LIB) flux.o gamma.o

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o
