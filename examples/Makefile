SHELL = /bin/sh

CC = gcc

INCLUDE_DIR = /usr/include
LIB_DIR = /usr/lib

BIN_PARAMS = params
BIN_TEST = test

CCFLAGS = -O2 -c -g
CLFLAGS = -o

%.o : %.c
	${CC} -I$(INCLUDE_DIR) $(CCFLAGS) $< -o $@

all: test.o params.o
	$(CC) params.c -I$(INCLUDE_DIR) -L$(LIB_DIR) -lm -lcparamlib $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -I$(INCLUDE_DIR) -L$(LIB_DIR) -lm -lcparamlib $(CLFLAGS) $(BIN_TEST)

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.o $(BIN_PARAMS) $(BIN_TEST)
