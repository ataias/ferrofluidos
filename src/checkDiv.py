#!/usr/bin/env python3

import glob, os, re, shutil
from vplot import *
from numpy import *
import time as tt
import math
from pathlib import Path

#Use this for a non-staggered grid
def divInside(Fx, Fy, n):
    divMax = 0.0
    dx = 1/n
    for i in range(2,n-1):
        for j in range(2,n-1):
            #Atenção! O índice j aqui é a direção x e o índice i a direção y
            #Quando o python leu do arquivo, esses índices foram trocados!
            divF = (Fx[i,j+1]-(Fx[i,j-1]))/(2*dx) + (Fy[i+1,j]-(Fy[i-1,j]))/(2*dx)
            if abs(divF) > abs(divMax):
                divMax = abs(divF)
    return divMax
    #end divInside

#Calcula rotacional nos pontos internos e retorna o maior valor de rotacional
def rotInside(Fx, Fy, n):
    rotMax = 0.0
    dx = 1/n
    for i in range(2,n-1):
        for j in range(2,n-1):
            rotF = (Fy[i,j+1]-Fy[i,j-1])/(2*dx) - (Fx[i+1,j]-Fx[i-1,j])/(2*dx)
            if abs(rotF) > abs(rotMax):
                rotMax = abs(rotF)

    return rotMax

def plotBEvolution(t, maxdivB, sideText, filename):
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

    plot(t, maxdivB)
    grid(True)
    axis([min(t), max(t), min(maxdivB)*1.1, max(maxdivB)*1.1])
    sideText()
    savefig(filename, dpi=200)

if __name__ == "__main__":

    p = Path('.')
    directories = [x for x in p.iterdir() if x.is_dir()]

    def Btext():
        xlabel(r'$t$')
        ylabel(r'max $|\nabla\cdot B|$')
        title(r'max $|\nabla\cdot B|$ evolution')

    def Htext():
        xlabel(r'$t$')
        ylabel(r'max $|\nabla\cdot H|$')
        title(r'max $|\nabla\cdot H|$ evolution')

    def Mtext():
        xlabel(r'$t$')
        ylabel(r'max $|\nabla\cdot M|$')
        title(r'max $|\nabla\cdot M|$ evolution')

    def Vtext():
        xlabel(r'$t$')
        ylabel(r'max $|\nabla\cdot \mathbf{v}|$')
        title(r'max $|\nabla\cdot \mathbf{v}|$ evolution')

    def Hrottext():
        xlabel(r'$t$')
        ylabel(r'max $|\nabla\times \mathbf{H}|$')
        title(r'max $|\nabla\times \mathbf{H}|$ evolution')


    for directoryPath in directories:
        directory = str(directoryPath)
        os.chdir(directory)

        parameters = re.findall(r"[-+]?\d*\.\d+|\d+", directory)

        Re = float(parameters[0])
        n = int(parameters[1])
        Pe = float(parameters[2])
        alpha = float(parameters[3])
        Cpm = float(parameters[4])
        t = float(parameters[5])
        fps = int(parameters[6])

        print("Pasta: ", directory)
        print("Re = ", Re)
        print("n = ", n)
        print("Pe = ", Pe)
        print("alpha = ", alpha)
        print("Cpm = ", Cpm)
        print("t = ", t)
        print("fps = ", fps)

        #Abrir arquivo
        f = open(directory + ".dat", 'rb')

        tvector = linspace(0, t, int(t*fps+1))
        maxdivV = zeros(size(tvector))
        maxdivB = zeros(size(tvector))
        maxdivH = zeros(size(tvector))
        maxdivM = zeros(size(tvector))
        maxrotH = zeros(size(tvector))

        #Como há o frame 0, precisa somar 1 ao número de frames
        numberFrames = int(round(fps*t))+1

        u = zeros((n,n), dtype=float64)
        v = zeros((n,n), dtype=float64)
        Hx = zeros((n,n), dtype=float64)
        Hy = zeros((n,n), dtype=float64)
        Mx = zeros((n,n), dtype=float64)
        My = zeros((n,n), dtype=float64)

        for i in range(0, numberFrames):
            time = i/fps
            # batchRead(u, v, p, Hx, Hy, Mx, My, phi, angles, n, f)
            #1 = SEEK_CUR -> offset relative to current position

            readMatrix(u, f, n)
            readMatrix(v, f, n)
            f.seek(n*n*8*1, 1)
            readMatrix(Hx, f, n)
            readMatrix(Hy, f, n)
            readMatrix(Mx, f, n)
            readMatrix(My, f, n)

            Bx = Hx+Mx
            By = Hy+My

            maxdivV[i] = divInside(u, v, n)
            maxdivB[i] = divInside(Bx, By, n)
            maxdivH[i] = divInside(Hx, Hy, n)
            maxdivM[i] = divInside(Mx, My, n)
            maxrotH[i] = rotInside(Hx, Hy, n)

            f.seek(n*n*8*2, 1)
        f.close()

        plotBEvolution(tvector,maxdivV, Vtext, "divV" + directory + ".png")
        plotBEvolution(tvector,maxdivB, Btext, "divB" + directory + ".png")
        plotBEvolution(tvector,maxdivH, Htext, "divH" + directory + ".png")
        plotBEvolution(tvector,maxdivM, Mtext, "divM" + directory + ".png")
        plotBEvolution(tvector,maxrotH, Hrottext, "rotH" + directory + ".png")


        os.chdir("..")
