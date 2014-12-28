#!/usr/bin/env python3

from numpy import *
from pylab import *
import struct
import sys

#Example of use: compare.py 41 1.0 1764.dat

#argv[1] is the size of matriz
#argv[2] is the time of simulation
#argv[3] is the fortran .dat file
#the name of the file from julia is obtained with argv[1]

n = int(sys.argv[1])
dx = 1/n
tot = int((n-1)**2)
uv = zeros((tot,8),dtype=float64)
#x, y, u, v (fortran), u, v (ataias), diff u, diff v

u = zeros((n-1,n-1), dtype=float64)
v = zeros((n-1,n-1), dtype=float64)
vort = zeros((n-1,n-1), dtype=float64)

#with open('1764.dat') as f:
with open(sys.argv[3]) as f:
    for line in f:
        numbers_str = line.split()
        #convert numbers to floats
        numbers_float = [float(x) for x in numbers_str]
        if numbers_float:
          i = round(numbers_float[0]/dx) - 1
          j = round(numbers_float[1]/dx) - 1
          vort[i,j] = numbers_float[3]
          u[i,j]    = numbers_float[4]
          v[i,j]    = numbers_float[5]

#use LaTeX, choose nice some looking fonts and tweak some settings
matplotlib.rc('font', family='serif')
matplotlib.rc('font', size=16)
matplotlib.rc('legend', fontsize=16)
matplotlib.rc('legend', numpoints=1)
matplotlib.rc('legend', handlelength=1)
matplotlib.rc('legend', frameon=False)
matplotlib.rc('xtick.major', pad=7)
matplotlib.rc('xtick.minor', pad=7)
matplotlib.rc('text', usetex=True)
matplotlib.rc('text.latex', 
             preamble=[r'\usepackage[T1]{fontenc}',
                       r'\usepackage{amsmath}',
                       r'\usepackage{txfonts}',
                       r'\usepackage{textcomp}'])

close('all')
figure(figsize=(8, 8))

x=linspace(0, 1 - 2*dx, n-1) + (dx)
y=linspace(0, 1 - 2*dx, n-1) + (dx)

x, y=meshgrid(x, y)
# plot vector field
vx = u.transpose()
vy = v.transpose()
Q = quiver(x, y, vx, vy, pivot='middle', headwidth=4, headlength=6)
title("Resultado - codigo da Camila")
# qk = quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
xlabel('$x$')
ylabel('$y$')
axis([0, 1.02, 0, 1.02])
show()

#Da forma abaixo posso ler os valores de um arquivo
# de pontos flutuantes
f = open('N' + sys.argv[1] + '.dat', 'rb')
t = float(sys.argv[2]) #time of simulation

numberFrames = round(180*t)
#f.seek((numberFrames - 1)*n*n*8*3) #get steady state solution
f.seek(-n*n*8*3, 2)
ua = zeros((n,n), dtype=float64)
va = zeros((n,n), dtype=float64)

for j in range(n):
	for i in range(n):
		ua[i,j] = struct.unpack('d',f.read(8))[0]

for j in range(n):
	for i in range(n):
		va[i,j] = struct.unpack('d',f.read(8))[0]

f.close()

um = zeros((n-1,n-1), dtype=float64)
vm = zeros((n-1,n-1), dtype=float64)
k = 0
for i in range(1,n):
  for j in range(1,n):
    um[i-1,j-1] = (ua[i-1,j-1] + ua[i,j-1] + ua[i-1,j] + ua[i,j])/4
    vm[i-1,j-1] = (va[i-1,j-1] + va[i,j-1] + va[i-1,j] + va[i,j])/4
    print("Ponto = ", k)
    k = k + 1

figure(figsize=(8, 8))
vxx = um.transpose()
vyy = vm.transpose()
title("Resultado - Ataias")
Q = quiver(x, y, vxx, vyy, pivot='middle', headwidth=4, headlength=6)
# qk = quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
xlabel('$x$')
ylabel('$y$')
axis([0, 1.02, 0, 1.02])
show()
# print(uv)

figure(figsize=(8, 8))
vxx = um.transpose()
vyy = vm.transpose()
title("Diferenca")
Q = quiver(x, y, vxx - vx, vyy, pivot='middle', headwidth=4, headlength=6)
# qk = quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
xlabel('$x$')
ylabel('$y$')
axis([0, 1.02, 0, 1.02])
show()

print("Maior diferença em u é ", abs(vxx - vx).max())
print("Maior valor em u - fortran é, ", abs(vx).max())
print("Maior valor em u - ataias é, ", abs(vxx).max())
print("Maior diferença em v é ", abs(vyy - vy).max())
print("Maior valor em v - fortran é, ", abs(vy).max())
print("Maior valor em v - ataias é, ", abs(vyy).max())
