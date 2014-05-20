#Autor: Ataias Pereira Reis
#Make file para projeto ferrofluidos
#Notice: need to create folders "bin" and "lib" before running this!
#Criado em: 19 de Maio de 2014
CC=clang++
EIGEN= ~/opt
PROJECT_HEADERS=./include
CXXFLAGS= -g -Wall -I $(EIGEN) -I $(PROJECT_HEADERS)
SRC= $(wildcard *.cc)
OBJ=$(SOURCES:.cc=.o)
VPATH=bin/lib:src:bin
vpath %.cc src/

#Nomes dos programas executáveis com funções main
# aparecem como dependências de "all"
all:	dir poisson navierstokes

#Seção com comandos para compilar cada programa

poisson:	poisson.cc sparsePD.cc
			$(CC) $(CXXFLAGS) $^ -o $@
			@mv $@ bin/$@  #move arquivos executáveis
			@mv $@.dSym bin/

navierstokes:	navierstokes.cc naviertest.cc
				$(CC) $(CXXFLAGS) $^ -o $@
#				@mv $? bin/lib #move dependências de arquivos modificados, outras já estão em lib
				@mv $@ bin/$@  #move arquivos executáveis
				@mv $@.dSym bin/

#Apaga arquivos objetos
clean:
	@cd bin/ && rm -rf *

#Apaga arquivos executáveis
mrproper:	clean

dir:
	@mkdir -p bin/