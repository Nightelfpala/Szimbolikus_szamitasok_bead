# Makefile

CC = g++

OUT = futtathato

all:
	$(CC) output.cpp -std=c++11 -o $(OUT) -L/usr/lib/ -lflint -lflint-arb -lmpfr

preproc:
	$(CC) test.cpp -std=c++11 -E > output.cpp