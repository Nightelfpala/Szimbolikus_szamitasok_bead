# Makefile

CC = g++

OUT = futtathato

all:
	$(CC) -o $(OUT) test.cpp -std=c++11 -L/usr/lib/ -lflint -lflint-arb -lmpfr