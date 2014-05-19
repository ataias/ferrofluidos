<a name="Vortex"/>
## Vortex
Este projeto faz parte do grupo Vortex.
- **Vortex:** <http://www.vortex.unb.br>

<a name="Ferrofluidos"/>
## Ferrofluidos

Este é um Projeto de Iniciação Científica na área de fluidos realizado na Universidade de Brasília sob a supervisão principal do professor Yuri Dumaresq. O objetivo é simular fluidos magnéticos em diversas aplicações. Para isso, discretizações de equações diferenciais devem ser realizadas e então executadas em um computador. As equações principais discretizadas aqui são: Laplace, Poisson e Navier Stokes (esta última ainda está em projeto). A solução das equações se dá num quadrado de lado 1 com o número de pontos na malha, n, variável. As discretizações são feitas aqui com o método das diferenças finitas.

<a name="Dependências"/>
## Dependências

Uma das dependências é a biblioteca Eigen. Outra é um compilador para C++, dois bastante utilizados são o GNU GCC e o LLVM Clang. Para criar gráficos, a biblioteca matplotlib para python será utilizada, mas isto ainda não foi implementado.

- **Eigen:** <http://eigen.tuxfamily.org>
- **GCC:** <https://gcc.gnu.org>
- **Clang:** <http://clang.llvm.org>

<a name="Como compilar e executar"/>
## Compilar e executar

Os programas aqui são feitos em C++ e os gráficos são feitos com a ajuda da biblioteca python matplotlib. 

Primeiramente, adquira o código fonte clonando o repositório git:

	git clone git@github.com:ataias/ferrofluidos.git

Após isso, entre no diretório `ferrofluidos/` e crie as pastas para os binários serem alocados, então execute `make`.

	cd ferrofluidos/
	mkdir bin && mkdir bin/lib
	make

Note que para este comando ser capaz de funcionar, tanto o compilador C++ como a biblioteca Eigen devem estar no caminho padrão para o sistema os poder encontrar. Em relação à Eigen, pode-se modificar arquivo `Makefile` no diretório `ferrofluidos` indicando onde ela está instalada, em qualquer seção do sistema. Por exemplo, obtenha uma versão da Eigen, como a 3.2.1:

	wget -nc http://bitbucket.org/eigen/eigen/get/3.2.1.tar.bz2
	tar xf 3.2.1.tar.bz2
	mkdir ~/opt #pode escolher outro diretório, e caso exista, não precisa desta linha
	mv eigen-eigen-6b38706d90a9/Eigen ~/opt/Eigen

Com a Eigen e outras bibliotecas localizadas no diretório `~/opt`, o compilador irá considerá-las além daquelas que estão no caminho padrão. 