module Magnetism

using Poisson
using NavierStokes
using NavierTypes

export getPhi!, getH!, getM!, getForce!, H0, solve!, Angle!

function R(theta)
    return [cos(theta) -sin(theta); sin(theta) cos(theta)]
end

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
end #end divInside

#takes values of Y vectors in the staggered grid to obtain an average
# in the same point of an X vector

function getAverageYtoX(i, j, A)
  return (A[i,j] + A[i,j+1] + A[i-1,j] + A[i-1,j+1])/4
end

#takes values of x vectors in the staggered grid to obtain an average
# in the same point of a Y vector
function getAverageXtoY(i,j, A)
  return (A[i,j] + A[i+1,j] + A[i+1,j-1] + A[i,j-1])/4
end

function getNextM!(M_next::VF, M::VF, M0::VF, H::VF, u, v, c1, c2, dt, dx)
  Mxt = 0.0; Myt = 0.0
  Hxt = 0.0; Hyt = 0.0

  for i in 2:n-1
    for j in 2:n-1
      Myt = getAverageYtoX(i,j, M.y)
      Hyt = getAverageYtoX(i,j, H.y)

      M_next.x[i,j]  = M.x[i,j] - c1*dt*(M.x[i,j]-M0.x[i,j])
      M_next.x[i,j] += -c2*dt*Myt*(M.x[i,j]*Hyt - Myt*H.x[i,j])
      M_next.x[i,j] += -0.25*dt*Myt*(v[i,j]+v[i,j+1]-v[i-1,j]-v[i-1,j+1])/dx
      M_next.x[i,j] += +0.25*dt*Myt*(u[i,j+1]-u[i,j-1])/dx
    end
  end #end for para M_next.x

  for i in 2:n-1
    for j in 2:n-1
      Mxt = getAverageXtoY(i,j, M.x)
      Hxt = getAverageXtoY(i,j, H.x)

      M_next.y[i,j] = M.y[i,j]
      M_next.y[i,j] -= c1*dt*(M.y[i,j]-M0.y[i,j])
      M_next.y[i,j] += c2*dt*Mxt*(Mxt*H.y[i,j] - M.y[i,j]*Hxt)
      M_next.y[i,j] += 0.25*dt*Mxt*(v[i+1,j]-v[i-1,j])/dx
      M_next.y[i,j] += 0.25*dt*Mxt*(-u[i,j]-u[i+1,j]+u[i,j-1]+u[i+1,j-1])/dx
    end
  end #end for para M_next.y

end #end getNextM()

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
end #end rotInside

# Calcula ângulo entre dois vetores para as matrizes H e M
#cada ponto deve ser um número complexo
#as entradas devem ser non-staggered
function Angle!(Hx, Hy, Mx, My, angles, n)
    for i in 1:n
        for j in 1:n
            angles[i,j] = angle(Hx[i,j],Hy[i,j], Mx[i,j], My[i,j])
        end
    end
    return angles;
end

function angle(hx, hy, mx, my)
  theta = atan2(hy,hx)
  x, y = R(theta)'*[mx; my]
  return atan2(y,x)*180/pi
end



function H0(Br, mu_m, L0, a, b, theta, x0, y0)
    K = (Br/pi/mu_m)
    H0x = x -> K*(atan(a*b/(x*sqrt(a^2+b^2+x^2))) - atan(a*b/((x+L0)*sqrt(a^2+b^2+(x+L0)^2))))

    Hx = (x,y) -> H0x(x*cos(theta)+y*sin(theta))*cos(theta)
    Hy = (x,y) -> H0x(x*cos(theta)+y*sin(theta))*sin(theta)
    return ((x,y) -> Hx(x-x0,y-y0)), ((x,y) -> Hy(x-x0,y-y0))
end # end H0()

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

    #apesar de ser calculado de 1 a n, somente pontos 2 a n-1 são usados
    # na rotina de Poisson
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

function getH!(n, phi, Hx, Hy)
    dx = 1/(n-2)
    for i in 2:n-1
        for j in 2:n-1
            Hx[i,j] = -(phi[i,j] - phi[i-1,j])/dx
            Hy[i,j] = -(phi[i,j] - phi[i,j-1])/dx
        end
    end

    #H fora da cavidade estava sendo calculado por meio de
    for j in 2:n-1
    #     Hy[1,j] = -(phi[1,j] - phi[1,j-1])/dx
    #     Hy[n,j] = -(phi[n,j] - phi[n,j-1])/dx
    # end
    #
    # for i in 2:n-1
    #     Hx[i,1] = -(phi[i,1] - phi[i-1,1])/dx
    #     Hx[i,n] = -(phi[i,n] - phi[i-1,n])/dx
    # end
    # Mas isto não é necessário. Como H é utilizado para obter M e M é nulo fora da cavidade
end #end getH!

function getM!(n, chi, phi, Mx, My, Hx, Hy)
  # Mx[:,:] = chi*Hx[:,:]
  # My[:,:] = chi*Hy[:,:]
  #Shliomis...
  x = 0 # variável auxiliar para calcular Mx0
  y = 0 # variável auxiliar para calcular My0

  for i in 2:n-1
    for j in 2:n-1

      x = sqrt(Hx[i,j])
    end
  end
end #end getM!


function getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy)
    dx = 1/(n-2)
    Myt = 0.0
    Mxt = 0.0
    Mb = 0.0
    Mc = 0.0

    #For internal points
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

    #For border points
    #Fronteira esquerda
    i = 2
    for j in 2:n-1
      #ao invés dessa média, que usa a malha fantasma
      # Myt = (My[i,j]+My[i,j+1]+My[i-1,j]+My[i-1,j+1])/4
      #usa-se uma interpolação com pontos internos
      # b é My no meio da primeira célula
      Mb = (My[i,j]+My[i,j+1])/2
      # c é My no meio célula à direta da primeira
      Mc = (My[i+1,j]+My[i+1,j+1])/2
      # Myt é uma interpolação linear a partir das últimas células mencionadas
      Myt = Mb + (Mc - Mb) * (-1/2)
      fx[i,j]  = Cpm*Mx[i,j]*( Hx[i,j+1]-Hx[i,j-1])/(2*dx)
      fx[i,j] += Cpm*Myt*(Hx[i+1,j]-Hx[i-1,j])/(2*dx)
    end

    #Fronteira direita
    i = n
    for j in 2:n-1
      #ao invés dessa média, que usa a malha fantasma
      # Myt = (My[i,j]+My[i,j+1]+My[i-1,j]+My[i-1,j+1])/4
      #usa-se uma interpolação com pontos internos
      # Mb é My no meio da célula atual
      Mb = (My[i-1,j]+My[i-1,j+1])/2
      # Mc é My no meio célula à direta da atual
      Mc = (My[i-2,j]+My[i-2,j+1])/2
      # Myt é uma interpolação linear para obter My no local de Mx
      # usando Mb e Mc
      Myt = Mb + (Mc - Mb) * (-1/2)
      fx[i,j]  = Cpm*Mx[i,j]*( Hx[i,j+1]-Hx[i,j-1])/(2*dx)
      fx[i,j] += Cpm*Myt*(Hx[i+1,j]-Hx[i-1,j])/(2*dx)
    end

    #Fronteira inferior
    j = 2
    for i in 2:n-1
      #usa-se uma interpolação com pontos internos
      # Mb é Mx no meio da célula atual
      Mb = (Mx[i,j]+My[i+1,j])/2
      # Mc é Mx no meio célula à direta da atual
      Mc = (Mx[i,j+1]+My[i,j+1])/2
      # Mxt é uma interpolação linear para obter Mx no local de My
      # usando Mb e Mc
      Mxt = Mb + (Mc - Mb) * (-1/2)
      fy[i,j]  = Cpm*Mxt*(Hy[i+1,j]-Hy[i-1,j])/(2*dx)
      fy[i,j] += Cpm*My[i,j]*( Hy[i,j+1]-Hy[i,j-1])/(2*dx)
    end

    #Fronteira superior
    j = n
    for i in 2:n-1
      #usa-se uma interpolação com pontos internos
      # Mb é Mx no meio da célula atual
      Mb = (Mx[i,j-1]+My[i+1,j-1])/2
      # Mc é Mx no meio célula à direta da atual
      Mc = (Mx[i,j-2]+My[i,j-2])/2
      # Mxt é uma interpolação linear para obter Mx no local de My
      # usando Mb e Mc
      Mxt = Mb + (Mc - Mb) * (-1/2)
      fy[i,j]  = Cpm*Mxt*(Hy[i+1,j]-Hy[i-1,j])/(2*dx)
      fy[i,j] += Cpm*My[i,j]*( Hy[i,j+1]-Hy[i,j-1])/(2*dx)
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

    hxn = zeros(n-2, n-2)
    hyn = zeros(n-2, n-2)
    mxn = zeros(n-2, n-2)
    myn = zeros(n-2, n-2)
    phin = zeros(n-2, n-2)
    angles = zeros(n-2,n-2)
    staggered2not!(Hx, Hy, phi, hxn, hyn, phin, n)
    staggered2not!(Mx, My, phi, mxn, myn, phin, n)
    Angle(hxn, hyn, mxn, myn, angles, n-2)
    file = open("problem0.dat", "w");
    write(file, mxn)
    write(file, myn)
    write(file, hxn)
    write(file, hyn)
    write(file, phin)
    write(file, angles)
    close(file)

    print("Max div B = ", divInside(mxn+hxn, myn+hyn, n-2), "\n")
    print("Max rot H = ", rotInside(hxn, hyn, n-2), "\n")
    Cpm = 10;
    fx = zeros(n,n);
    fy = zeros(n,n);
    getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy);

end

end
