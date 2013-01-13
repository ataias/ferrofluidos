# coding: utf-8
from dolfin import *
# Create mesh and define function space
mesh = UnitSquare(6, 4)
V = FunctionSpace(mesh, "Lagrange", 1)
# Define boundary conditions
u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
def u0_boundary(x, on_boundary):
    return 0 #sem condições de dirichlet
bc = DirichletBC(V, u0, u0_boundary)
#Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
a = inner(nabla_grad(u), nabla_grad(v))*dx
g = Constant(0)
L = f*v*dx - g*v*ds #g é zero... não é grande coisa
# Compute solution
u = Function(V)
solve(a == L, u, bc)
# Plot solution and mesh
plot(u)
plot(mesh)
# Dump solution to file in VTK format
file = File("poisson.pvd")
file << u
# Hold plot
interactive()
