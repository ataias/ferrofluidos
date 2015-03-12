module Magnetism

using Poisson
using NavierStokes

export getPhi!, getMH!, getForce!, H0

function H0(Br, mu_m, L0, a, b, theta, x0, y0)
    K = (Br/pi/mu_m)
    H0x = x -> K*(atan(a*b/(x*sqrt(a^2+b^2+x^2))) - atan(a*b/((x+L0)*sqrt(a^2+b^2+(x+L0)^2))))
    
    Hx = (x,y) -> H0x(x*cos(theta)+y*sin(theta))*cos(theta)
    Hy = (x,y) -> H0x(x*cos(theta)+y*sin(theta))*sin(theta)
    return ((x,y) -> Hx(x-x0,y-y0)), ((x,y) -> Hy(x-x0,y-y0))
end

#fHx and fHy are functions
function getPhi!(n, phi, Mx, My, fHx, fHy, A)
    dx = 1/(n-2)
    
    #f = div M
    f = zeros(n,n); #valor de f no meio da célula
    for i in 2:n-1
        for j in 2:n-1
            f[i,j] = (Mx[i+1,j] - Mx[i,j])/dx + (My[i,j+1] - My[i,j])/dx
        end
    end
    
    left = zeros(n);
    right = zeros(n);
    upper = zeros(n);
    lower = zeros(n);
    
    for i in -1:n-2
      x = (i+0.5)*dx;
      left[i+2]  = Mx[2,i+2] - fHx(0,x)
      right[i+2] = Mx[n,i+2] - fHx(1,x)
      upper[i+2] = My[i+2,n] - fHy(x,1)
      lower[i+2] = My[i+2,2] - fHy(x,0)
    end
    #Solve for p with Neumann conditions here
    poissonNeumannSparseSolver(n, phi, f, A, left, right, upper, lower)
    
end #end getPhi!

function getMH!(n, chi, phi, Mx, My, Hx, Hy)
    dx = 1/(n-2)
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
    dx = 1/(n-2)
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
    dx = 1/(n-2);
    phi = zeros(n,n);
    Mx = zeros(n,n);
    My = zeros(n,n);
#    
    fMx = (x,y) -> x
    fMy = (x,y) -> -y
    
    for i in -1:n-2
        x = (i+0.5)*dx;
        for j in -1:n-2
            y = (j+0.5)*dx;
            Mx[i+2,j+2] = fMx(x,y)
            My[i+2,j+2] = fMy(x,y)
        end
    end
    
    Br = 1.2;
    mu_0 = 4pi*1e-7;
    mu_m = mu_0;
#    L0 = 0.0254
    L0 = 2
    a, b = 0.025, 0.025;
    x0, y0 = -0.2, 0;
    theta = 0;
#    (fHx, fHy) = H0(Br, mu_m, L0, a, b, theta, x0, y0);
    A = getANeumannSparse(n);
    getPhi!(n, phi, Mx, My, (x,y) -> 1, (x,y) -> 1, A);
    
    Hx = zeros(n,n);
    Hy = zeros(n,n);
    getMH!(n, chi, phi, Mx, My, Hx, Hy);
    
    Cpm = 10;
    fx = zeros(n,n);
    fy = zeros(n,n);
    getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy);
    
end

end