#!/usr/bin/env python3

from matplotlib import *
from pylab import *
import struct
import sys

#Modules to list files of current module and remove filename extensions as needed
import glob, os, re

def readMatrix(array, file, n):
    for i in range(n):
    	for j in range(n):
            array[i,j] = struct.unpack('d',file.read(8))[0]

#file "f"
def batchRead(u, v, p, Hx, Hy, Mx, My, phi, angles, n, f):
    readMatrix(u, f, n)
    readMatrix(v, f, n)
    readMatrix(p, f, n)
    readMatrix(Hx, f, n)
    readMatrix(Hy, f, n)
    readMatrix(Mx, f, n)
    readMatrix(My, f, n)
    readMatrix(phi, f, n)
    readMatrix(angles, f, n)

def plotStreamFrame(u, v, x, y, n, sideText, time, filename):
    #use LaTeX, choose nice some looking fonts and tweak some settings
    rc('font', family='serif')
    rc('font', size=16)
    rc('legend', fontsize=16)
    rc('legend', numpoints=1)
    rc('legend', handlelength=1)
    rc('legend', frameon=False)
    rc('xtick.major', pad=7)
    rc('xtick.minor', pad=7)
    rc('text', usetex=True)
    rc('text.latex',
                 preamble=[r'\usepackage[T1]{fontenc}',
                           r'\usepackage{amsmath}',
                           r'\usepackage{txfonts}',
                           r'\usepackage{textcomp}'])

    close('all')
    figure(figsize=(12, 8))

    #1 é a velocidade máxima na malha, ela varia se a condição de contorno for modificada
    norm = Normalize(vmin=0, vmax=1.0) #intensidade deve ser de 0 a 1
    streamplot(x, y, u, v, color=u, linewidth=1.3, cmap=cm.brg, arrowsize=4, norm=norm)
    colorbar(norm=norm, cmap=cm.winter, ticks=[0, 0.25, 0.75, 1])
    xlabel('$x$')
    ylabel('$y$')
    sideText(time)

    axis([0, 1.0, 0, 1.02])
    savefig(filename, dpi=200)

def plotPressure(x, y, p, filename):
    #Plotar pressão
    close('all')
    fig = figure()
    maxl = max(p)
    minl = min(p)
    levels = arange(minl, maxl, (maxl - minl)/10.0)
    CS = contour(x, y, p)
    clabel(CS, inline=1, fontsize=10)
    title('Pressure')
    savefig(filename, dpi=200)

def plotAngle(x, y, angles, n, filename):
  close('all')
  fig = figure()
  dx = dy = 1.0/n #quadrado de lado 1
  cmap = get_cmap('PiYG')
#  levels = MaxNLocator(nbins=15).tick_values(angles.min(), angles.max())
  levels = MaxNLocator(nbins=30).tick_values(-180.0, 180.0)
  contourf(x, y, angles, levels=levels, cmap=cmap)
  colorbar()
  title(r'Angle between $\mathbf{H}_{\mathrm{calc}}$ and $\mathbf{M}$')
  xlabel(r'$x$')
  ylabel(r'$y$')
  savefig(filename, dpi=200)

def plotPointEvolution(t, w, sideText, filename):
    #use LaTeX, choose nice some looking fonts and tweak some settings
    rc('font', family='serif')
    rc('font', size=16)
    rc('legend', fontsize=16)
    rc('legend', numpoints=1)
    rc('legend', handlelength=1)
    rc('legend', frameon=False)
    rc('xtick.major', pad=7)
    rc('xtick.minor', pad=7)
    rc('text', usetex=True)
    rc('text.latex',
                 preamble=[r'\usepackage[T1]{fontenc}',
                           r'\usepackage{amsmath}',
                           r'\usepackage{txfonts}',
                           r'\usepackage{textcomp}'])

    close('all')
    figure(figsize=(12, 8))

    plot(t, w)

    sideText()

    axis([min(t), max(t), min(w), max(w)*1.1])
    grid(True)
    savefig(filename, dpi=200)

def plotMEvolution(t, modM, phaseM, phaseDiffMH, sideText, filename):
    #use LaTeX, choose nice some looking fonts and tweak some settings
    rc('font', family='serif')
    rc('font', size=16)
    rc('legend', fontsize=16)
    rc('legend', numpoints=1)
    rc('legend', handlelength=1)
    rc('legend', frameon=False)
    rc('xtick.major', pad=7)
    rc('xtick.minor', pad=7)
    rc('text', usetex=True)
    rc('text.latex',
                 preamble=[r'\usepackage[T1]{fontenc}',
                           r'\usepackage{amsmath}',
                           r'\usepackage{txfonts}',
                           r'\usepackage{textcomp}'])

    close('all')
    figure(figsize=(12, 8))

    subplot(2,1,1)
    plot(t, modM)
    grid(True)
    axis([min(t), max(t), min(modM), max(modM)*1.1])
    ylabel('$|\mathbf{M}|$')
    sideText()


    subplot(2,1,2)
    plot(t,phaseM, label='ang(M)')
    plot(t,phaseDiffMH, label='ang(H-M)')
    xlabel('$t$')
    ylabel('Phase (degrees)')

    grid(True)
    legend()
    savefig(filename, dpi=200)


if __name__ == "__main__":

    # makePNGforVideo('N' + str(n) + '.dat', n)
    # f = open('N' + str(n) + '.dat', 'rb')
    # t = float64(sys.argv[2])
    # step = int(sys.argv[3])
    # chi = float64(sys.argv[4])
    # Cpm = float64(sys.argv[5])
    # Re = float64(sys.argv[6])
    # gamma = float64(sys.argv[7])

    n = 100
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

    batchRead(u, v, p, Hx, Hy, Mx, My, phi, angles, n, f)

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
