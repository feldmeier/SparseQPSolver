
PATH_LAPACK = ../LAPACK

all: playground

playground:
	gcc -c denseNASM.c -o playground.o -std=c99
	gcc -c parser.c -o parser.o -std=c99
	gcc -o denseNASM playground.o parser.o $(PATH_LAPACK)/liblapack.a $(PATH_LAPACK)/liblapacke.a $(PATH_LAPACK)/librefblas.a -lgfortran -lm
	rm playground.o
clean:
	rm *.o 


