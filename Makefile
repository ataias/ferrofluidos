#Autor: Ataias Pereira Reis
#Make file para projeto ferrofluidos
#Notice: need to create folders "bin" and "lib" before running this!
#Criado em: 19 de Maio de 2014
CC=clang++
EIGEN= ~/opt
PROJECT_HEADERS=./include
CXXFLAGS=-c -Wall -I $(EIGEN) -I $(PROJECT_HEADERS)
SRC= $(wildcard *.cc)
OBJ=$(SOURCES:.cc=.o)
VPATH=bin/lib:src:bin

#Nomes dos programas executáveis com funções main
# aparecem como dependências de "all"
all:	poisson

#Seção que faz o linking cada programa de acordo com
#os objetos que ele necessita
#move os objetos pra pasta bin/lib
poisson:	poisson.o sparsePD.o
			$(CC) $^ -o $@
			@mv $? bin/lib #move dependências de arquivos modificados, outras já estão em lib
			@mv $@ bin/$@  #move arquivos executáveis

#Seção que compila arquivos, mas sem fazer o link
$(OBJ):	$(SRC)
	$(CC) $(CXXFLAGS) $(OBJ)

#Apaga arquivos objetos
clean:
	@cd bin/lib && rm -rf *.o
	@cd bin/ && rm -rf *.o

#Apaga arquivos executáveis
mrproper:	clean
	@cd bin && rm -rf poisson