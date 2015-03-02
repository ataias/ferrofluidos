module Magnetism

using Poisson
using NavierStokes

function getPhi!(n, phi, Mx, My, left, right, upper, lower)
    dx = 1/(n-2)
    
    #f = div M
    f = zeros(n,n); #valor de f no meio da célula
    for i in 2:n-1
        for j in 2:n-1
            f[i,j] = (Mx[i+1,j] - Mx[i,j])/dx + (My[i,j+1] - My[i,j])/dx
        end
    end
    
    solveDirichletPoissonExplicit!(n, phi, f, left, right, upper, lower)
end #end getPhi!

function getMH!(n, chi, phi, Mx, My, Hx, Hy)
    for i in 2:n-1
        for j in 2:n-1
            Hx[i,j] = -(phi[i,j] - phi[i-1,j])/dx
            Hy[i,j] = -(phi[i,j] - phi[i,j-1])/dx
        end
    end
    
    for j in 3:n-1
        Hy[1,j] = -(phi[1,j] - phi[1,j-1])/dx
        Hy[n,j] = -(phi[n,j] - phi[n,j-1])/dx
    end
    
    for i in 3:n-1
        Hx[i,1] = -(phi[i,1] - phi[i-1,1])/dx
        Hx[i,n] = -(phi[i,n] - phi[i-1,n])/dx
    end
    
    Mx[:,:] = chi*Hx[:,:]
    My[:,:] = chi*Hy[:,:]
end #end getMH!


function getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy)
    Myt = 0.0
    Mxt = 0.0
    
    for i in 2:n-1
        for j in 2:n-1
            Myt = (My[i,j]+My[i,j+1]+My[i-1,j]+My[i-1,j+1])/4
            fx[i,j]  = Cpm*Mx[i,j]*( Hx[i,j+1]-Hx[i,j-1])/(2*dx)
            fx[i,j] += Cpm*Myt*(Hx[i+1,j]-Hx[i-1,j])/(2*dx)
            
            Mxt = (Mx[i,j]+Mx[i+1,j]+Mx[i,j-1]+Mx[i+1,j-1])/4
            fy[i,j]  = Cpm*Mxt*(Hy[i+1,j]-Hy[i-1,j])/(2*dx)
            fy[i,j] += Cpm*My[i,j]*( Hy[i,j+1]-Hy[i,j-1])/(2*dx)
        end
    end 

end #end getForce!

function solve!(n = 7)
    dx = 1/(n-2)
    phi = zeros(n,n);
    Mx = zeros(n,n);
    My = zeros(n,n);
    
    left = zeros(1,n);
    right = zeros(1,n);
    upper = zeros(1,n);
    lower = zeros(1,n);
    for i in -1:n-2
      x = (i+0.5)*dx;
      upper[i+2] = (sinpi(x))^2;
    end
    
    getPhi!(n, phi, Mx, My, left, right, upper, lower)
    
    chi = 2
    Hx = zeros(n,n)
    Hy = zeros(n,n)
    getMH!(n, chi, phi, Mx, My, Hx, Hy)
    
    Cpm = 4
    fx = zeros(n,n);
    fy = zeros(n,n);
    
end

end