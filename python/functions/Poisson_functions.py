#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys
sys.path.append("/home/ataias/git/ff/build/src/") #localização das bibliotecas
sys.path.append("/home/ataias/git/ff/python/functions/") #localização das bibliotecas
from Poisson_functions import *
from libpoisson import Poisson #biblioteca com os códigos para resolver a equação de Poisson
import scipy
import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D #Para plotar gráficos 3D
import matplotlib.pyplot as plt
from numpy import linalg as LA #linear algebra, norm function
import time

def plot_dSolution(SIZE, dSolution):
    pylab.figure(1)
    ax = pylab.subplot(111, projection='3d')
    
    x = pylab.linspace(0,1,SIZE)
    xx, yy = pylab.meshgrid(x,x)
    zz = dSolution[0:SIZE,0:SIZE]
    
    ax.plot_surface(xx, yy, zz,
                    rstride=1,
                    cstride=1,
                    cmap=pylab.cm.coolwarm,
                    linewidth=0.2)
    
    ax.set_xlabel(r'$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    pylab.title(r'$\nabla^2u = 1$ Dirichlet conditions')
    #pylab.savefig('plot3d.pdf')
    
    pylab.figure(2)
    ax1 = pylab.subplot(111)
    
    x = pylab.linspace(0,1,SIZE)
    xx1, yy1 = pylab.meshgrid(x,x)
    #zz = yy**2/2
    #zz = dSolution[1:SIZE-1,1:SIZE-1]
    cs = ax1.contour(xx1,yy1,zz,15)
    ax1.grid(True)
    ax1.axis('scaled')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel('$y$')
    pylab.clabel(cs,inline=1,fontsize=10)
    pylab.title(r'Level curves of the function')
    #pylab.savefig('plot3d.png')
    pylab.show()

def calculateError():
    dError = np.array([[0 for j in range(0,SIZE)] for i in range(0,SIZE)],dtype=np.float64)
    dDeltaX = 1.0/(SIZE-1)
    dDeltaX2 = 1.0/(SIZE-1)**2
    for i in range(0,SIZE):
        for j in range(0,SIZE):
            #QUINAS
            if i==0 and j==0: 
                continue;
            if i==0 and j==SIZE-1: 
                continue;
            if i==SIZE-1 and j==0:
                continue;
            if i==SIZE-1 and j==SIZE-1:
                continue;
        
            if i==0:        #NORTH POINT
                continue;
            elif i==SIZE-1: #SOUTH POINT
                continue;
            elif j==0:      #WEST POINT
                continue;
            elif j==SIZE-1: #EAST POINT
                continue;
            else:
                dError[i,j] = (dSolution[i-1,j]+dSolution[i+1,j]+\
                dSolution[i,j-1]+dSolution[i,j+1]-4*dSolution[i,j])/dDeltaX2 - dNonHomogeneity[i,j]
    return LA.norm(dError)

def solvePoisson(SIZE, dSolution, cMeshPoisson):
    #Condições de Fronteira
    cMesh = cMeshPoisson;
    dBoundaryConditions = np.array([[0 for j in range(0,SIZE)] for i in range(0,SIZE)],dtype=np.float64)
    dBoundaryConditions[0,:]=[0 for j in range(0,SIZE)];
    dBoundaryConditions[:,0]=[75 for j in range(0,SIZE)];
    dBoundaryConditions[SIZE-1,:]=[50 for j in range(0,SIZE)];
    dBoundaryConditions[:,SIZE-1]=[0 for j in range(0,SIZE)];
    #Não homogeneidade
    dNonHomogeneity = np.array([[0 for j in range(0,SIZE)] for i in range(0,SIZE)],dtype=np.float64)
    #Especificando o tipo de problema
    bDirichletOrNeumann = 0 #Indica que o método utilizado será de Dirichlet nas quatro fronteiras
    bSparseOrNot = 0 #Indica que o método será direto
    nMatrixOrder = SIZE #Indica ordem da matriz
    #Inicializando programa em C++
    cMesh.conditions(dBoundaryConditions,
                     dNonHomogeneity,
                     bDirichletOrNeumann,
                     bSparseOrNot,
                     nMatrixOrder);
    #Obtendo solução
    cMesh.poissonSolver();
    cMesh.saveSolution(dSolution);
