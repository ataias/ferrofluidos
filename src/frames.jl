using Transient
using NavierTypes

#modo de usar
#n   - ARGS[1] é o tamanho da matriz escalonada
#t   - ARGS[2] é o tempo de simulação, em segundos
#Re  - ARGS[3] é o número de Reynolds
#divFactor  - ARGS[4] dividir o passo de tempo, deve ser maior do que 1
#chi - ARGS[5]
#Cpm - ARGS[6]
#save - ARGS[7] se for 1, salva um arquivo com as matrizes evoluindo no tempo
#a - ARGS[8] é o deslocamento em x da posição central do campo magnético
#b - ARGS[9] é o deslocamento em y da posição central do campo magnético
#gamma - ARGS[10] é a intensidade do campo magnético
#Exemplo:
# 		julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 -0.05 3
# A saída padrão é salva num arquivo txt nomeado de acordo com Re e n

n = int(ARGS[1]);
t = float(ARGS[2]);
Re = float(ARGS[3]);
divFactor = float(ARGS[4]);
dt = getDt(n, Re, float(ARGS[4]));
chi = float(ARGS[5]);
Cpm = float(ARGS[6]);
save = bool(int(ARGS[7]));
a = float(ARGS[8]);
b = float(ARGS[9]);
gamma = float(ARGS[10]);

filename = "Re" * ARGS[3] * "N" * string(n-2) *".txt"
file = open(filename, "w")
redirect_stdout(file)
@time transient(n, dt, Re, t, chi, Cpm, gamma, a, b, divFactor, save)
close(file)