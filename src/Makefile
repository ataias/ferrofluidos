#Autor: Ataias Pereira Reis
#Make file para projeto ferrofluidos
CC=clang++
EIGEN= ~/opt
PROJECT_HEADERS=../include
CFLAGS=-c -Wall -I $(EIGEN) -I $(PROJECT_HEADERS)
#LDFLAGS=
#SOURCES=poisson.cc sparsePD.cc
#OBJECTS=$(SOURCES:.cc=.o)
#VPATH = src
#BUILDDIR = build
#EXECUTABLE=poisson

all:	poisson

poisson:	poisson.o sparsePD.o
			$(CC) poisson.o sparsePD.o -o poisson

poisson.o:	poisson.cc
	$(CC) $(CFLAGS) poisson.cc -o poisson.o

sparsePD.o: sparsePD.cc
			$(CC) $(CFLAGS) sparsePD.cc -o sparsePD.o

clean:
	rm -rf *.o

mrproper:	clean
	rm -rf poisson