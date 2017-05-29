# Makefile

CC = g++

OUT = futtathato

all:
	$(CC) -o $(OUT) test.cpp -std=c++11 -L/usr/lib/ -lflint -lflint-arb -lmpfr
	
arbpp_test:
	$(CC) -o $(OUT) test_arbpp.cpp -std=c++11 -L/usr/lib/ -lflint -lflint-arb -lmpfr