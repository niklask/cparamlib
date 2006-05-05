SHELL = /bin/sh

CC = gcc
AR = ar
SWIG = /opt/bin/swig

LIB = libcparammodel.a
BIN_PARAMS = params
BIN_TEST = test

PYTHON_LIB = _cparammodel.so
PYTHON_INCLUDE = /opt/include/python2.4

CCFLAGS = -c -g
CFLAGS  = -c
CLFLAGS = -o
ARFLAGS = rcs
SWIGFLAGS = -python -o

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

py: lib
	$(SWIG) $(SWIGFLAGS) cparammodel_wrap.c cparammodel.i
	$(CC) -shared cparammodel_wrap.c -L. -lm -lcparammodel -I$(PYTHON_INCLUDE) $(CLFLAGS) $(PYTHON_LIB)

lib: test.o sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o params.o
	$(AR) $(ARFLAGS) $(LIB) sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o
	$(CC) params.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_TEST)

all: py

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o *.py* $(BIN_PARAMS) $(BIN_TEST) $(PYTHON_LIB)
