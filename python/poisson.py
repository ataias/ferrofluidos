# coding: utf-8
#!/usr/bin/env python

"""
Biblioteca que contém classes e métodos para resolver a equação de Poisson 
com condições de Dirichlet ou Neumann.

Autor: Ataias Pereira Reis <ataiasreis at gmail.com>
Última modificação: Março, 2013
"""

import numpy
from scipy import weave
from scipy.weave import converters
import sys, time

#Global variables
DIRICHLET = 'Dirichlet';
METHOD = 'Gauss'
    
class Grid:
    
    """Uma classe de malha para salvar detalhes básicos."""
    
    def __init__(self, nx=31, ny=31, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.u = numpy.zeros((nx, ny), 'd')
        # Utilizado para calcular erro
        self.old_u = self.u.copy()        

    def setBC(self, l, r, b, t):        
        """Configura a condição de contorno dados os valores da esquerda, direita, de baixo
        e de cima (ou arrays)"""        
        self.u[0, :] = l
        self.u[-1, :] = r
        self.u[:, 0] = b
        self.u[:,-1] = t
        self.old_u = self.u.copy()

    def setBCFunc(self, func):
        """Configura a BC de acordo com uma função de duas variáveis dada."""
        xmin, ymin = self.xmin, self.ymin
        xmax, ymax = self.xmax, self.ymax
        x = numpy.arange(xmin, xmax + self.dx*0.5, self.dx)
        y = numpy.arange(ymin, ymax + self.dy*0.5, self.dy)
        self.u[0 ,:] = func(xmin,y)
        self.u[-1,:] = func(xmax,y)
        self.u[:, 0] = func(x,ymin)
        self.u[:,-1] = func(x,ymax)

    def computeError(self):        
        """Computes absolute error using an L2 norm for the solution.
        This requires that self.u and self.old_u must be appropriately
        setup."""        
        v = (self.u - self.old_u).flat
        return numpy.sqrt(numpy.dot(v,v))
    

class PoissonSolver:
    
    """Um resolvePoisson.""" 
    global DIRICHLET, METHOD
    def __init__(self, 
                grid,
                method = METHOD,
                cLeft = DIRICHLET, 
                cRight = DIRICHLET, 
                cBottom= DIRICHLET,
                cTop = DIRICHLET
                ):
        self.grid = grid;
        self.setSolver(method);

    def inlineGauss(self, dt=0.0):        
        """Takes a time step using inlined C code -- this version uses
        blitz arrays."""        
        g = self.grid
        nx, ny = g.u.shape
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u
        
        code = """
               double tmp, err, diff;
               err = 0.0;
               for (int i=1; i<nx-1; ++i) {
                   for (int j=1; j<ny-1; ++j) {
                       tmp = u(i,j);
                       u(i,j) = ((u(i-1,j) + u(i+1,j))*dy2 +
                                 (u(i,j-1) + u(i,j+1))*dx2)*dnr_inv;
                       diff = u(i,j) - tmp;
                       err += diff*diff;
                   }
               }
               return_val = sqrt(err);
               """
        # compiler keyword only needed on windows with MSVC installed
        err = weave.inline(code,
                           ['u', 'dx2', 'dy2', 'dnr_inv', 'nx','ny'],
                           type_converters = converters.blitz,
                           compiler = 'gcc')
        return err         
    
    def setSolver(self, method='Gauss'):
        if method == 'Gauss':
            self.solverMethod = self.inlineGauss;
            
    
    def solve(self, n_iter=0, eps=1.0e-16):        
        """Solves the equation given an error precision -- eps.  If
        n_iter=0 the solving is stopped only on the eps condition.  If
        n_iter is finite then solution stops in that many iterations
        or when the error is less than eps whichever is earlier.
        Returns the error if the loop breaks on the n_iter condition
        and returns the iterations if the loop breaks on the error
        condition."""        
        err = self.solverMethod()
        count = 1

        while err > eps:
            if n_iter and count >= n_iter:
                return err
            err = self.solverMethod()
            count = count + 1

        return count


def BC(x, y):    
    """Used to set the boundary condition for the grid of points.
    Change this as you feel fit."""    
    return (x**2 - y**2)


def test(nx=5, ny=5, eps=1.0e-16, n_iter=0, method='Gauss'):
    g = Grid(nx, ny)
    g.setBCFunc(BC)
    s = PoissonSolver(g, method)
    t = time.clock()
    s.solve(n_iter=n_iter, eps=eps)
    print s.grid.u.transpose()
    return time.clock() - t
    

def main(n=31, n_iter = 0):
    test();

if __name__ == "__main__":
    main()