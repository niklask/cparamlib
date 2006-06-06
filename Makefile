SHELL = /bin/sh

CC = gcc
AR = ar
SWIG = /opt/bin/swig

LIB = libcparammodel.a
BIN_PARAMS = params
BIN_TEST = test

PYTHON_LIB = _cparammodel.so
PYTHON_PATH = python/
PYTHON_INCLUDE = /opt/include/python2.4

CCFLAGS = -O2 -c -g
CFLAGS  = -c
CLFLAGS = -o
ARFLAGS = rcs
SWIGFLAGS = -python

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

py: lib
	$(SWIG) $(SWIGFLAGS) -outdir $(PYTHON_PATH) cparammodel.i
	$(CC) -shared cparammodel_wrap.c -L. -lm -lcparammodel -I$(PYTHON_INCLUDE) $(CLFLAGS) $(PYTHON_PATH)$(PYTHON_LIB)

lib: test.o sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o params.o
	$(AR) $(ARFLAGS) $(LIB) sigma.o gamma.o elec.o posi.o nue.o numu.o antinue.o antinumu.o
	$(CC) params.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_TEST)

all: py

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.a *.o *_wrap.c $(BIN_PARAMS) $(BIN_TEST) $(PYTHON_PATH)$(PYTHON_LIB) $(PYTHON_PATH)cparammodel.py
