SHELL = /bin/sh

CC = gcc
AR = ar

LIB = libcparammodel.a

CCFLAGS = -c -g
CLFLAGS = -o
ARFLAGS = rcs

all: cparamlib

%.o : %.cpp
	${CC} $(CCFLAGS) $< -o $@

cparamlib: flux.o gamma.o
	$(AR) $(ARFLAGS) $(LIB) flux.o gamma.o

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o
