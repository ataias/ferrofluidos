using Transient
using NavierTypes

#modo de usar
#n   - ARGS[1] é o tamanho da matriz escalonada
#t   - ARGS[2] é o tempo de simulação, em segundos
#Re  - ARGS[3] é o número de Reynolds
#divFactor  - ARGS[4] dividir o passo de tempo, deve ser maior do que 1
#c1 - ARGS[5]
#Cpm - ARGS[6]
#save - ARGS[7] se for 1, salva um arquivo com as matrizes evoluindo no tempo
#a - ARGS[8] é o deslocamento em x da posição central do campo magnético
#b - ARGS[9] é o deslocamento em y da posição central do campo magnético
# alpha - ARGS[10] é relacionado com a intensidade do campo magnético
#Exemplo:
# 		julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 -0.05 3
# (nomag)     julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 -0.05 3
# A saída padrão é salva num arquivo txt nomeado de acordo com Re e n


#TODO
#Algo que me incomoda aqui é que há algumas coisas que podem levar a erros
# e que não são verificadas
#divFactor deve ser maior que um, mas não uso nenhum "assert" pra ver isso,
# ou um if mesmo
# mesmo coisa com o tempo de simulação, deve ser positivo
#a e b devem ser números que indiquem posições fora da fronteira. A singularidade
# passando por dentro do escoamento é algo que não funciona direito

#o python não pede para dizer o 2 em "52" por exemplo, posso editar isso aqui também
#fiz algumas modificações aí

n = int(ARGS[1]);
@assert(n > 10, "Arg 1: Use um tamanho de malha maior.");

t = float(ARGS[2]);
if t <= 0
   println("Arg 2: Tempo de simulação deve ser maior que zero.")
   exit()
end

Re = float(ARGS[3]);
if Re <= 0
  println("Arg 3: Re deve ser maior que zero.")
  exit()
end

divFactor = float(ARGS[4]);
if divFactor <= 1
  println("Arg 4: Fator de divisão deve ser maior que 1.")
  exit()
end

dt = getDt(n, Re, float(ARGS[4]));

c1 = float(ARGS[5]);
Cpm = float(ARGS[6]);

save = bool(int(ARGS[7]));

a = float(ARGS[8]);
b = float(ARGS[9]);
if (0 < b < 1) && (0 < a < 1)
  println("Arg 8 'a' e 9 'b': Centro do campo magnético deve ser definido fora da cavidade.")
  exit()
end

alpha = float(ARGS[10]);


filename = "Re" * string(int(Re)) * "N" * string(n-2) *".txt"
file = open(filename, "w")
redirect_stdout(file)
@time transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1)
close(file)
