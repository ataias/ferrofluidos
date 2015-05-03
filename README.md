<a name="Vortex"/>
## Vortex
Este projeto faz parte do grupo Vortex.
- **Vortex:** <http://www.vortex.unb.br>

<a name="Ferrofluidos"/>
## Ferrofluidos

Este é um Projeto de Iniciação Científica na área de fluidos realizado na Universidade de Brasília sob a supervisão principal do professor Yuri Dumaresq. O objetivo é simular fluidos magnéticos em diversas aplicações. Para isso, discretizações de equações diferenciais devem ser realizadas e então executadas em um computador. As equações principais discretizadas aqui são: Laplace, Poisson e Navier Stokes (esta última ainda está em projeto). A solução das equações se dá num quadrado de lado 1 com o número de pontos na malha, n, variável. As discretizações são feitas aqui com o método das diferenças finitas.

<a name="Dependências"/>
## Dependências

- **Julia:** <http://julialang.org>
- **Python 3:** <https://www.python.org>
- **Matploblib:** <http://matplotlib.org>
- **NumPy:** <http://www.numpy.org>
- **LaTeX:** <http://latex-project.org>

<a name="Como compilar e executar"/>
## Compilar e executar

Os programas aqui são feitos em Julia e os gráficos são feitos com a ajuda da biblioteca python matplotlib. Para o código em Julia, basta executar e o sistema fará a compilação Just-In-Time. No caso do python, os códigos serão interpretados.

Primeiramente, obtenha o código fonte clonando o repositório git. Isto pode ser feito tanto por ssh como por https:

	git clone git@github.com:ataias/ferrofluidos.git
	git clone https://github.com/ataias/ferrofluidos.git

Após isso, entre no diretório `ferrofluidos/` e então terás acesso ao código fonte. Julia usa compilação Just-In-Time. Para executar o código `frames.jl`, pode-se ir ao terminal e:

	julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 0.0 > out.txt

Detalhes sobre os parâmetros estão disponíveis no código `frames.jl`. Dados da velocidade e pressão estão disponíveis na saída `.dat` de cada um desses programas. Pode-se processar a saída do `frames.jl` com o `vectorField.py`

	chmod +x vectorField.py
	./vectorField.py

É recomendável adicionar a linha 

	push!(LOAD_PATH,"/home/ataias/ferrofluidos/src")

ao arquivo `~/.juliarc.jl`, lembrando de adicionar o diretório correto onde `ferrofluidos/src` está.
