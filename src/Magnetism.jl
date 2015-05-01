module Magnetism

using Poisson
using NavierStokes

export getPhi!, getMH!, getForce!, H0, solve!, Angle!

#Use this for a non-staggered grid
function divInside(Fx, Fy, n)
    divMax = 0.0
    dx = 1/n
    for i in 2:n-1
        for j in 2:n-1
            divF = (Fx[i+1,j]-(Fx[i-1,j]))/(2*dx) + (Fy[i,j+1]-(Fy[i,j-1]))/(2*dx)
            if abs(divF) > abs(divMax)
                divMax = divF
            end
        end
    end
    return divMax
end

#use this for a non-staggered grid
function rotInside(Fx, Fy, n)
    dx = 1/n
    rotMax = 0.0
    for i in 2:n-1
        for j in 2:n-1
            rotF = (Fy[i+1,j]-(Fy[i-1,j]))/(2*dx) - (Fx[i,j+1]-(Fx[i,j-1]))/(2*dx)
            if abs(rotF) > abs(rotMax)
                rotMax = rotF
#                print("rot(", i, ",", j, ") =", rotF, "\n")
            end
        end
    end
    return rotMax
end

# Calcula ângulo entre dois vetores para as matrizes H e M
#cada ponto deve ser um número complexo
function Angle(a, b)
    angles = zeros(size(a,1), size(a,2));
    for i in 1:size(a, 1)
        for j in 1:size(a,2)
            angles[i,j] = angle(reim(a[i,j]),reim(b[i,j]))
        end
    end
    return angles;
end

function angle(a, b)
  return acosd(dot(a,b)/norm(a)/norm(b))
end

function angle(a::(Number,Number), b::(Number,Number))
    return angle([a[1], a[2]], [b[1], b[2]])
end


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

    for j in 2:n-1
        Hy[1,j] = -(phi[1,j] - phi[1,j-1])/dx
        Hy[n,j] = -(phi[n,j] - phi[n,j-1])/dx
    end

    for i in 2:n-1
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

    theta = pi/4
    R = theta -> [cos(theta) -sin(theta); sin(theta) cos(theta)]
    Rt = R(theta)

    for i in -1:n-2
        x = (i+0.5)*dx;
        for j in -1:n-2
            y = (j+0.5)*dx;
            xl, yl = Rt'*[x y]'
            Mx[i+2,j+2] = (Rt*[fMx(xl,yl) fMy(xl,yl)]')[1]
            My[i+2,j+2] = (Rt*[fMx(xl,yl) fMy(xl,yl)]')[2]
        end
    end
    #Quando for calcular phi, devo fazer M = chi*H primeiro? Isso terá influência no phi por meio de suas condições de contorno.

    A = getANeumannSparse(n);
    getPhi!(n, phi, Mx, My, (x,y) -> 1, (x,y) -> 1, A);

    chi = 1
    Hx = zeros(n,n);
    Hy = zeros(n,n);
    getMH!(n, chi, phi, Mx, My, Hx, Hy);

    hxn = zeros(n-2, n-2); hyn = zeros(n-2, n-2);
    mxn = zeros(n-2, n-2); myn = zeros(n-2, n-2);
    phin = zeros(n-2, n-2);
    staggered2not!(Hx, Hy, phi, hxn, hyn, phin, n);
    staggered2not!(Mx, My, phi, mxn, myn, phin, n);
    file = open("problem0.dat", "w");
    write(file, mxn)
    write(file, myn)
    write(file, hxn)
    write(file, hyn)
    write(file, phin)
    close(file)

    print("Max div B = ", divInside(mxn+hxn, myn+hyn, n-2), "\n")
    print("Max rot H = ", rotInside(hxn, hyn, n-2), "\n")
    Cpm = 10;
    fx = zeros(n,n);
    fy = zeros(n,n);
    getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy);

end

end
