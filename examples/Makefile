SHELL = /bin/sh

CC = gcc

BIN_PARAMS = params
BIN_TEST = test

CCFLAGS = -O2 -c -g
CLFLAGS = -o

%.o : %.c
	${CC} $(CCFLAGS) $< -o $@

all: test.o params.o
	$(CC) params.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_PARAMS)
	$(CC) test.c -L. -lm -lcparammodel $(CLFLAGS) $(BIN_TEST)

clean:
	-@ echo "Cleaning up..."
	-@ rm -r *.o $(BIN_PARAMS) $(BIN_TEST)
