SHELL = /bin/sh

CC = gcc

INCLUDE_DIR = ${HOME}/usr/include
LIB_DIR = ${HOME}/usr/lib

BIN_PARAMS = params
BIN_TEST = test

CCFLAGS = -O2 -I$(HOME)/usr/include -c -g
CLFLAGS = -I$(HOME)/usr/include -o

%.o : %.c
	${CC} -I$(INCLUDE_DIR) $(CCFLAGS) $< -o $@

all: test.o params.o
	$(CC) params.c -I$(INCLUDE_DIR) -L$(LIB_DIR) -lm -lcparamlib $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -I$(INCLUDE_DIR) -L$(LIB_DIR) -lm -lcparamlib $(CLFLAGS) $(BIN_TEST)

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.o $(BIN_PARAMS) $(BIN_TEST)
