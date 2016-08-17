#!/usr/bin/env python3

from matplotlib import *
use('Agg')
from pylab import *
import struct
import sys

#Modules to list files of current module and remove filename extensions as needed
import glob, os, re

def setup():
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

def plotStreamFrame(u, v, x, y, n, sideText, time, filename):
    setup()

    #1 é a velocidade máxima na malha, ela varia se a condição de contorno for modificada
    norm = Normalize(vmin=0, vmax=1.0) #intensidade deve ser de 0 a 1
    streamplot(x, y, u, v, color=u, linewidth=1.3, cmap=cm.brg, arrowsize=4, norm=norm)
    colorbar(norm=norm, cmap=cm.winter, ticks=[0, 0.25, 0.75, 1])
    xlabel('$x$')
    ylabel('$y$')
    sideText(time)

    axis([0, 1.0, 0, 1.02])
    savefig(filename, format="png", dpi=200)

def plotVectorFrame(u, v, x, y, n, sideText, time, filename):
    setup()

    #1 é a velocidade máxima na malha, ela varia se a condição de contorno for modificada
    # norm = Normalize(vmin=0, vmax=1.0) #intensidade deve ser de 0 a 1
    # streamplot(x, y, u, v, color=u, linewidth=1.3, cmap=cm.brg, arrowsize=4, norm=norm)
    # colorbar(norm=norm, cmap=cm.winter, ticks=[0, 0.25, 0.75, 1])
    scale = (amax(u)**2 + amax(v)**2)**0.5
    quiver(x[0:n:4], y[0:n:4], u[0:n:4], v[0:n:4])
    xlabel('$x$')
    ylabel('$y$')
    sideText(time)

    axis([0, 1.0, 0, 1.02])
    savefig(filename, format="png", dpi=200)

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
    savefig(filename, format="png", dpi=200)

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
  savefig(filename, format="png", dpi=200)

def plotPointEvolution(t, w, sideText, filename):
    #use LaTeX, choose nice some looking fonts and tweak some settings
    setup()

    plot(t, w)

    sideText()

    axis([min(t), max(t), min(w), max(w)*1.1])
    grid(True)
    savefig(filename, format="png", dpi=200)

def plotMEvolution(t, modM, phaseM, phaseDiffMH, sideText, filename):
    setup()

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
