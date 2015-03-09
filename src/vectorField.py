#!/usr/bin/env python3
#Example to run this program:
#./vectorField.py 100 3.0 4 1.0 4 30
#arg[1] - mesh size
#arg[2] - time of simulation
#arg[3] - step (if size is n, then n/step vectors will be plotted, this makes things neat)
#arg[4] - chi
#arg[5] - Cpm
#arg[6] - Re    

# import useful modules
import matplotlib 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import math
import sys

from pylab import *
from numpy import *
import struct

def readFromFile(array, file, n):
    for i in range(n):
    	for j in range(n):
            array[i,j] = struct.unpack('d',file.read(8))[0]

#file "f"
def batchRead(u, v, p, Hx, Hy, phi, n, f):
    readFromFile(u, f, n)
    readFromFile(v, f, n)
    readFromFile(p, f, n)
    readFromFile(Hx, f, n)
    readFromFile(Hy, f, n)
    readFromFile(phi, f, n)

def plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, filename):
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
    Q = quiver(x[0:n:step,0:n:step], y[0:n:step,0:n:step], u[0:n:step,0:n:step], v[0:n:step,0:n:step], pivot='middle', headwidth=4, headlength=6)
    qk = quiverkey(Q, 0.5, 1.0, 1, r'$\mathbf{v}$, mesh $' + str(n) + r'\times' + str(n) + '$, $\chi = ' + str(chi) + '$, Cpm = ' + str(Cpm) + ', Re = ' + str(Re), fontproperties={'weight': 'bold'})
    xlabel('$x$')
    ylabel('$y$')
    axis([0, 1.0, 0, 1.02])
    savefig(filename, dpi=300)

def plotPressure(x, y, p, filename):
    #Plotar pressão
    close('all')
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(x, y, p, rstride=1, cstride=1, cmap=cm.coolwarm,
           linewidth=0, antialiased=False)
    pmin = int(amin(p)) - 1
    pmax = int(amax(p)) + 1
    ax.set_zlim(pmin, pmax)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig(filename, dpi=300)
    
#Calcula rotacional nos pontos internos e retorna o maior valor de rotacional
def rotInside(Fx, Fy, n):
    rotMax = 0.0
    for i in range(2,n-2):
        for j in range(2,n-2):
            rotF = (Fy[i+1,j]-(Fy[i-1,i]))/(2*dx) - (Fx[i,j+1]-(Fx[i,j-1]))/(2*dx)
            if abs(rotF) > abs(rotMax):
                rotMax = rotF
    return rotMax

def divInside(Fx, Fy, n):
    divMax = 0.0
    for i in range(2,n-2):
        for j in range(2,n-2):
            divF = (Fx[i+1,j]-(Fx[i-1,i]))/(2*dx) + (Fy[i,j+1]-(Fy[i,j-1]))/(2*dx)
            if abs(divF) > abs(divMax):
                divMax = divF
    return divMax

#To create video:
#ffmpeg -i png/%4d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p out.mp4
def makePNGforVideo(filename, n):
    f = open('N' + str(n) + '.dat', 'rb')
    t = float64(sys.argv[2])
    step = int(sys.argv[3])
    chi = float64(sys.argv[4])
    Cpm = float64(sys.argv[5])
    Re = float64(sys.argv[6])
    dx = 1/n
    
    u = zeros((n,n), dtype=float64)
    v = zeros((n,n), dtype=float64)
    p = zeros((n,n), dtype=float64)
    Hx = zeros((n,n), dtype=float64)
    Hy = zeros((n,n), dtype=float64)
    phi = zeros((n,n), dtype=float64)
    
    numberFrames = int(round(180*t))
    
    x=linspace(0, 1 - dx, n) + (dx/2)
    y=linspace(0, 1 - dx, n) + (dx/2)
    x, y=meshgrid(x, y)
    
    for i in range(0, numberFrames-2):
        batchRead(u, v, p, Hx, Hy, phi, n, f)
        plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, str(i).zfill(4) + ".png")
    
    f.close()
    
    
if __name__ == "__main__":
    
    n = int(sys.argv[1]) #esta é a dimensão da malha escalonada menos 2
    f = open('N' + str(n) + '.dat', 'rb')
#    makePNGforVideo('N' + str(n) + '.dat', n)
    t = float64(sys.argv[2])
    step = int(sys.argv[3])
    chi = float64(sys.argv[4])
    Cpm = float64(sys.argv[5])
    Re = float64(sys.argv[6])
    
    dx = 1/n
    numberFrames = round(180*t)
    #f.seek((numberFrames - 1)*n*n*8*3) #get steady state solution
    f.seek(-n*n*8*6, 2)
    u = zeros((n,n), dtype=float64)
    v = zeros((n,n), dtype=float64)
    p = zeros((n,n), dtype=float64)
    Hx = zeros((n,n), dtype=float64)
    Hy = zeros((n,n), dtype=float64)
    phi = zeros((n,n), dtype=float64)

    batchRead(u, v, p, Hx, Hy, phi, n, f)

    f.close()
    # generate grid
    x=linspace(0, 1 - dx, n) + (dx/2)
    y=linspace(0, 1 - dx, n) + (dx/2)
    x, y=meshgrid(x, y)
    plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, 'vectorField.png')
    plotVectorField(Hx, Hy, x, y, n, step, chi, Cpm, Re, 'vectorFieldH.png')
    
    plotPressure(x, y, p, 'pressure.png')

    ij = math.ceil((n+2)/2-1) - 1
    w = (v[ij,ij+1]-v[ij,ij-1])/(2*dx) - (u[ij+1,ij]-u[ij-1,ij])/(2*dx)
    print("Max rot inside is {0:.3f}k".format(rotInside(Hx, Hy, n)))
    print("Max div inside is {0:.3f}k".format(divInside(Hx, Hy, n)))
    #wc = -0.63925
    #error = abs(abs(wc)-abs(w))
    #print("Error in vorticity: ", error)

    #Da forma abaixo posso ler os valores de um arquivo
    # de pontos flutuantes
    #sys.argv[1] é o n
    

