#!/usr/bin/env python3

# import useful modules
import matplotlib 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from numpy import *
from pylab import *
import struct

#Da forma abaixo posso ler os valores de um arquivo
# de pontos flutuantes
f = open('N43.dat', 'rb')
t = 2.2
n = 43 #esta é a dimensão da malha escalonada menos 2
dx = 1/n

numberFrames = round(180*t)
#f.seek((numberFrames - 1)*n*n*8*3) #get steady state solu   tion
f.seek(-n*n*8*3, 2)
u = zeros((n,n), dtype=float64)
v = zeros((n,n), dtype=float64)
p = zeros((n,n), dtype=float64)

for i in range(n):
	for j in range(n):
		u[i,j] = struct.unpack('d',f.read(8))[0]

for i in range(n):
	for j in range(n):
		v[i,j] = struct.unpack('d',f.read(8))[0]

for i in range(n):
	for j in range(n):
		p[i,j] = struct.unpack('d',f.read(8))[0]

f.close()

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
 
# generate grid
x=linspace(0, 1 - dx, n) + (dx/2)
y=linspace(0, 1 - dx, n) + (dx/2)
x, y=meshgrid(x, y)
# calculate vector field
vx=u
vy=v
# plot vecor field
Q = quiver(x, y, vx, vy, pivot='middle', headwidth=4, headlength=6)
# qk = quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
xlabel('$x$')
ylabel('$y$')
axis([0, 1.0, 0, 1.02])
show()
# savefig('visualization_vector_fields_1.png')

#Plotar pressão
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, p, rstride=1, cstride=1, cmap=cm.coolwarm,
       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()



