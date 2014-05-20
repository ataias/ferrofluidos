#Autor: Ataias Pereira Reis
#Make file para projeto ferrofluidos
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
#executando make tenta criar o diretório bin/
#se ele já não existir
all:	dir poisson navierstokes

#Seção com comandos para compilar cada programa

poisson:	poisson.cc sparsePD.cc
			$(CC) $(CXXFLAGS) $^ -o $@
			@mv $@ bin/$@  #move arquivo executávei
			@mkdir -p $@.dSYM && mv $@.dSYM bin/

navierstokes:	navierstokes.cc naviertest.cc
				$(CC) $(CXXFLAGS) $^ -o $@
				@mv $@ bin/$@
				@mkdir -p $@.dSYM && mv $@.dSYM bin/

#Apaga arquivos em bin
clean:
	@cd bin/ && rm -rf *

#Apaga arquivos executáveis
mrproper:	clean

dir:
	@mkdir -p bin/