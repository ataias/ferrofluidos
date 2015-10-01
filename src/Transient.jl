module Transient

using NavierStokes
using NavierTypes
using Magnetism
using Poisson

export transient

##TODO
#Save Mx and My and process it correctly in .py

function factor(i)
    if i <= 500
        return i/500.0
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
function transient(n, dt, Re, t, Cpm, chi, a, b, save, fps, filename)
  dx = 1/(n-2)

  println("Dados sobre simulação: ")
  println(" n\t= ", n)
  println(" dx\t= ", dx)
  println(" t\t= ", t)
  println(" Re\t= ", Re)
  println(" dt\t= ", dt)
  println(" Cpm\t= ", Cpm)
  println(" χ\t= ", chi)
  println(" a\t= ", a)
  println(" b\t= ", b)
  println(strftime(time()), "\n")


	steps = integer(t/dt)
	c = integer(n/2); #center, non-staggered grid

  NS = createNSObject(n, Re)
	for i in -1:n-2
		NS.uB[i+2] = 0*(sinpi(i*dx))^2 #inicialmente, repouso
	end

  #Condições de contorno, basta editar em v_old
	NS.v_old.x[:,n] = 2*NS.uB

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
    #Salva valores das matrizes em t = 0
    if(save)
        write(file, un); write(file, vn);
        write(file, pn);
        write(file, Hxn); write(file, Hyn);
        write(file, Mxn); write(file, Myn);
        write(file, phin);
        write(file, angles);
    end

	for i in 1:steps
    fact = factor(i)
    for j in -1:n-2
        NS.uB[j+2] = fact*(sinpi(j*dx))^2
    end
    #Caso superparamagnético, é necessário calcular phi primeiro
    getPhi!(n, phi, Mx, My, (x,y) -> fact * fHx(x,y), (x,y) -> fact * fHy(x,y), A)
    getMH!(n, chi, phi, Mx, My, Hx, Hy)
    getForce!(n, Cpm, Hx, Hy, Mx, My, NS.f.x, NS.f.y);
    solve_navier_stokes!(NS)

		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(NS.v.x, NS.v.y, NS.p, un,  vn,  pn,   n)
            staggered2not!(Hx,     Hy,     phi,  Hxn, Hyn, phin, n)
            staggered2not!(Mx,     My,           Mxn, Myn,       n)
            Angle!(Hxn, Hyn, Mxn, Myn, angles, n-2)
            if(save)
                write(file, un); write(file, vn);
                write(file, pn);
                write(file, Hxn); write(file, Hyn);
                write(file, Mxn); write(file, Myn);
                write(file, phin);
                write(file, angles);
            end
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
			end
			F = simpson(tau, n-2)
			println("  F\t= ", F)
		end #if (i % timeToSave == 0) || (i == steps)

    #Preparing for next time step
		NS.v.x, NS.v_old.x = NS.v_old.x, NS.v.x
		NS.v.y, NS.v_old.y = NS.v_old.y, NS.v.y
	end

    if(save)
        close(file)
    end
	return 0
end

end
