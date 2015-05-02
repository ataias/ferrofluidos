#!/usr/bin/env python3
#Example to run this program:
#./vectorField.py 100 3.0 4 1.0 4 30
#arg[1] - mesh size
#arg[2] - time of simulation
#arg[3] - step (if size is n, then n/step vectors will be plotted, this makes things neat)
#arg[4] - chi
#arg[5] - Cpm
#arg[6] - Re    
#comando para fazer vídeo
#   convert -quality 100 *.png movie.mpg

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

#to plot angles
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator

def readFromFile(array, file, n):
    for i in range(n):
    	for j in range(n):
            array[i,j] = struct.unpack('d',file.read(8))[0]

#file "f"
def batchRead(u, v, p, Hx, Hy, phi, angles, n, f):
    readFromFile(u, f, n)
    readFromFile(v, f, n)
    readFromFile(p, f, n)
    readFromFile(Hx, f, n)
    readFromFile(Hy, f, n)
    readFromFile(phi, f, n)
    readFromFile(angles, f, n)

def plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, gamma, time, filename):
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
    figure(figsize=(12, 8))
    #Q = quiver(x[0:n:step,0:n:step], y[0:n:step,0:n:step], u[0:n:step,0:n:step], v[0:n:step,0:n:step], pivot='middle', headwidth=4, headlength=6)
    #qk = quiverkey(Q, 0.5, 1.0, 1, r'$\mathbf{v}$, mesh $' + str(n) + r'\times' + str(n) + '$, $\chi = ' + str(chi) + '$, Cpm = ' + str(Cpm) + ', Re = ' + str(Re), fontproperties={'weight': 'bold'})

    step = 1 #for streamplot, let step = 1 always
    streamplot(x[0:n:step,0:n:step], y[0:n:step,0:n:step], u[0:n:step,0:n:step], v[0:n:step,0:n:step], color=u, linewidth=1.3, cmap=cm.winter, arrowsize=4)
    colorbar()
    xlabel('$x$')
    ylabel('$y$')
    text(1.21, 1.0, r'Re = ' + str(Re))
    text(1.21, 0.95, r'$\chi = ' + str(chi) + '$')
    text(1.21, 0.9, r'Cpm = ' + str(Cpm))
    text(1.21, 0.85, r'$n = ' + str(n) + '$')
    text(1.21, 0.8, r'$\gamma = ' + str(gamma) + '$')
    title(r'$\mathbf{v}$, ' + 't = ' + '{:10.8f}'.format(time))

    axis([0, 1.0, 0, 1.02])
    savefig(filename, dpi=300)

def plotPressure(x, y, p, filename):
    #Plotar pressão
    close('all')
    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    surf = ax.plot_surface(x, y, p, rstride=1, cstride=1, cmap=cm.coolwarm,
#           linewidth=0, antialiased=False)
#    pmin = int(amin(p)) - 1
#    pmax = int(amax(p)) + 1
#    ax.set_zlim(pmin, pmax)
#    ax.zaxis.set_major_locator(LinearLocator(10))
#    ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
#    fig.colorbar(surf, shrink=0.5, aspect=5)
    CS = plt.contour(x, y, p)
    plt.clabel(CS, inline=1, fontsize=10)   
    plt.title('Pressure')
    plt.savefig(filename, dpi=300)
    
def plotAngle(x, y, angles, n, filename):
  close('all')
  fig = plt.figure()
  dx = dy = 1.0/n #quadrado de lado 1
  cmap = plt.get_cmap('PiYG')
#  levels = MaxNLocator(nbins=15).tick_values(angles.min(), angles.max())
  levels = MaxNLocator(nbins=15).tick_values(-180.0, 180.0)
  plt.contourf(x, y, angles, levels=levels, cmap=cmap)
  plt.colorbar()
  plt.title(r'Angle between $\mathbf{H}_{\mathrm{calc}}$ and $\mathbf{M}$')
  xlabel(r'$x$')
  ylabel(r'$y$')  
  plt.savefig(filename, dpi=300)
  
    
#Calcula rotacional nos pontos internos e retorna o maior valor de rotacional
def rotInside(Fx, Fy, n):
    rotMax = 0.0
    pos = (0, 0)
    dx = 1/n
    for i in range(2,n-2):
        for j in range(2,n-2):
            rotF = (Fy[i,j+1]-Fy[i,j-1])/(2*dx) - (Fx[i+1,j]-Fx[i-1,j])/(2*dx)
            if abs(rotF) > abs(rotMax):
                rotMax = rotF
                pos = (i, j)
                
    return rotMax, pos

def divInside(Fx, Fy, n):
    divMax = 0.0
    dx = 1/n
    for i in range(1,n-1):
        for j in range(1,n-1):
            divF = (Fx[i+1,j]-(Fx[i-1,j]))/(2*dx) + (Fy[i,j+1]-(Fy[i,j-1]))/(2*dx)
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
    gamma = float64(sys.argv[7])
    dx = 1/n
    
    u = zeros((n,n), dtype=float64)
    v = zeros((n,n), dtype=float64)
    p = zeros((n,n), dtype=float64)
    Hx = zeros((n,n), dtype=float64)
    Hy = zeros((n,n), dtype=float64)
    phi = zeros((n,n), dtype=float64)
    angles = zeros((n,n), dtype=float64)
    
    numberFrames = int(round(180*t))
    
    x=linspace(0, 1 - dx, n) + (dx/2)
    y=linspace(0, 1 - dx, n) + (dx/2)
    x, y=meshgrid(x, y)
    time = 0.0
    
    for i in range(0, numberFrames):
        time = (i)/180.0
        try:
            batchRead(u, v, p, Hx, Hy, phi, angles, n, f)
            plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, gamma, time, str(i).zfill(4) + ".png")
        except:
            break
    
    f.close()
    
    
if __name__ == "__main__":
    
    n = int(sys.argv[1]) #esta é a dimensão da malha escalonada menos 2
#    makePNGforVideo('N' + str(n) + '.dat', n)
    f = open('N' + str(n) + '.dat', 'rb')
    t = float64(sys.argv[2])
    step = int(sys.argv[3])
    chi = float64(sys.argv[4])
    Cpm = float64(sys.argv[5])
    Re = float64(sys.argv[6])
    gamma = float64(sys.argv[7])
    
    dx = 1/n
    numberFrames = round(180*t)
    f.seek(-n*n*8*7, 2)
    u = zeros((n,n), dtype=float64)
    v = zeros((n,n), dtype=float64)
    p = zeros((n,n), dtype=float64)
    Hx = zeros((n,n), dtype=float64)
    Hy = zeros((n,n), dtype=float64)
    phi = zeros((n,n), dtype=float64)
    angles = zeros((n,n), dtype=float64)

    batchRead(u, v, p, Hx, Hy, phi, angles, n, f)

    f.close()
    # generate grid
    x=linspace(0, 1 - dx, n) + (dx/2)
    y=linspace(0, 1 - dx, n) + (dx/2)
    
    x, y=meshgrid(x, y)
    
    plotPressure(x, y, p, 'pressure.png')
    plotAngle(x,y,angles,n, 'angle.png')
    
    plotVectorField(u, v, x, y, n, step, chi, Cpm, Re, gamma, t, 'vectorField.png')
    if(gamma>1e-15):
      plotVectorField(Hx, Hy, x, y, n, step, chi, Cpm, Re, gamma, t, 'vectorFieldH.png')
    

