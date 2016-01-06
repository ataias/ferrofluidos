module Transient

using NavierStokes
using NavierTypes
using Magnetism
using Poisson

export transient

##TODO
#Save Mx and My and process it correctly in .py
#I think this have been done, check

function factor(i)
    if i <= 100
        return i/100.0
    else
        return 1.0
    end
end

#transient
#simula t tempos adimensionais
#retorna o valor em cada tempo para o 0.5, 0.5
#resolve equações para um dado n e Re
#t is time, in dimensioless units, of physical simulation
#fps is frames per second
#convec é igual a 1.0 ou 0.0 e é multiplicado pelo termo convectivo para suprimí-lo
# ou permitir que seja considerado
function transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, filename, convec=0.0)
    dx = 1/(n-2)

    println("Dados sobre simulação: ")
    println(" n\t= ", n)
    println(" dx\t= ", dx)
    println(" t\t= ", t)
    println(" Re\t= ", Re)
    println(" dt\t= ", dt)
    println(" Cpm\t= ", Cpm)
    println(" α\t= ", alpha)
    println(" a\t= ", a)
    println(" b\t= ", b)
    println(" c1\t= ", c1)
    println(strftime(time()), "\n")


	steps = integer(t/dt)
	c = integer(n/2); #center, non-staggered grid

    NS = createNSObject(n, Re)

	for i in -1:n-2
		NS.uB[i+2] = 0*(sinpi(i*dx))^2 #inicialmente, repouso
	end

    #Condições de contorno, basta editar em v_old
	NS.v_old.x[:,n] = 2*NS.uB

    #Força
	#fx = zeros(n,n)
	#fy = zeros(n,n)

    if(save)
	   file = open(filename, "w")
    end

	un = zeros(n-2,n-2)+1e-15 #malha não escalonada
	vn = zeros(n-2,n-2)+1e-15
	pn = zeros(n-2,n-2)+1e-15

	numberFrames = integer(fps*t)
	timeToSave = integer(steps/numberFrames)

    #Variáveis para a parte magnética
    phi = zeros(n,n);
    Mx = zeros(n,n); Mx_old = zeros(n,n);
    My = zeros(n,n); My_old = zeros(n,n);
    left = zeros(n);
    right = zeros(n);
    upper = zeros(n);
    lower = zeros(n);
    for i in -1:n-2
      x = (i+0.5)*dx;
#      upper[i+2] = 1
#      left[i+2] = x
#      right[i+2] = x
       upper[i+2] = sinpi(x)^2
    end
    Hx = zeros(n,n)
    Hy = zeros(n,n)

    #non-staggered forms
    Hxn = zeros(n-2, n-2)+1e-15
    Hyn = zeros(n-2, n-2)+1e-15
    Mxn = zeros(n-2, n-2)+1e-15
    Myn = zeros(n-2, n-2)+1e-15
    phin = zeros(n-2, n-2)+1e-15
    A = getANeumannSparse(n);

    #gamma foi acoplado dentro do alpha que é usado para calcular M0
    fHx = (x,y) ->  1/(2*pi)*(y-b)/((x-a)^2+(y-b)^2)
    fHy = (x,y) -> -1/(2*pi)*(x-a)/((x-a)^2+(y-b)^2)

    fact = 0
    angles = zeros(n-2, n-2)+1e-15;
    fxn = zeros(n-2, n-2)+1e-15
    fyn = zeros(n-2, n-2)+1e-15
    
    #Salva valores das matrizes em t = 0
    if(save)
        write(file, un); write(file, vn);
        write(file, pn);
        write(file, Hxn); write(file, Hyn);
        write(file, Mxn); write(file, Myn);
        write(file, phin);
        write(file, angles);
        write(file, fxn)
        write(file, fyn)
    end

    #Magnetização em regime é constante e só depende de H aplicado
    #pode ser calculada só uma vez
    Mx0 = zeros(n,n)
    for i in 2:n #
      for j in 2:n-1
        x = (i-2)*dx
        y = (j-2)*dx + dx/2
        mH = sqrt(fHx(x,y)^2 + fHy(x,y)^2) #módulo de H
        Mx0[i,j] = (coth(alpha*mH) - 1/(alpha*mH)) * fHx(x,y)/mH
      end
    end

    My0 = zeros(n,n)
    for i in 2:n-1
      for j in 2:n
        x = (i-2)*dx + dx/2
        y = (j-2)*dx
        mH = sqrt(fHx(x,y)^2 + fHy(x,y)^2) #módulo de H
        My0[i,j] = (coth(alpha*mH) - 1/(alpha*mH)) * fHy(x,y)/mH
      end
    end

  #Variáveis auxiliares
  #Termo de convecção magnética
  v∇Mx = zeros(n,n)
  v∇My = zeros(n,n)

	for i in 1:steps
    fact = factor(i)
    for j in -1:n-2
        NS.uB[j+2] = fact*(sinpi(j*dx))^2
    end
    #O primeiro passo é obter M em cada passo de tempo
    #faz sentido não usar o fator fact para getM, pois a evolução inicia
    #do valor anterior de M, que é 0 no tempo 0
    v∇M!(n, NS.v.x, NS.v.y, Mx_old, My_old, v∇Mx, v∇My)
    getM!(n, c1, dt, Mx, My, Mx_old, My_old, Mx0, My0, convec*v∇Mx, convec*v∇My)
    getPhi!(n, phi, Mx, My, fHx, fHy, A)
    getH!(n, phi, Hx, Hy)
    getForce!(n, Cpm, Hx, Hy, Mx, My, NS.f.x, NS.f.y);
    solve_navier_stokes!(NS)

    #println("i = ", i, " and timeToSave = ", timeToSave, ", therefore i % timeToSave = ", i % timeToSave == 0)
		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(NS.v.x, NS.v.y, NS.p, un,  vn,  pn,   n)
      staggered2not!(Hx,     Hy,     phi,  Hxn, Hyn, phin, n)
      staggered2not!(Mx,     My,           Mxn, Myn,       n)
      staggered2not!(NS.f.x, NS.f.y,       fxn, fyn,       n)
      Angle!(Hxn, Hyn, Mxn, Myn, angles, n-2)
      if(save)
        write(file, un); write(file, vn);
        write(file, pn);
        write(file, Hxn); write(file, Hyn);
        write(file, Mxn); write(file, Myn);
        write(file, phin);
        write(file, angles);
        write(file, fxn);
        write(file, fyn);
      end #if(save)

			println("t = ", i*dt)
			println("  u[0.5,0.5]\t= ", un[c,c])
			println("  v[0.5,0.5]\t= ", vn[c,c])
      #Observe que vn e un estão salvas nas quinas.
      #Este cálculo de vorticidade só é válido caso n seja PAR!
			vortc =  ((vn[c+1,c]-vn[c-1,c]) - (un[c,c+1]-un[c,c-1]) )/(2*dx)
			println("  ω [0.5,0.5]\t= ", vortc)
            println("  Pressure values in range\t [", minimum(pn), ", ", maximum(pn), "]")
			tau = zeros(n-2)
			for i in 2:n-1
				tau[i-1] = (1/Re)*((NS.v.x[i,n]-NS.v.x[i,n-1])/dx)
			end #for i in 2:n-1

			F = simpson(tau, n-2)
			println("  F\t= ", F)
		end #if (i % timeToSave == 0) || (i == steps)

    #Preparing for next time step
		NS.v.x, NS.v_old.x = NS.v_old.x, NS.v.x
		NS.v.y, NS.v_old.y = NS.v_old.y, NS.v.y

    Mx, Mx_old = Mx_old, Mx
    My, My_old = My_old, My

	end #for i in 1:steps

  if(save)
    close(file)
  end
	return 0
end #function transient

end #module
