module Transient

using HDF5
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

macro saveFrames(index, dt)
  return quote
    frames["velocity/x/" * string($index)] = un
    frames["velocity/y/" * string($index)] = vn
    frames["pressure/" * string($index)] = pn
    frames["H/x/" * string($index)] = Hxn
    frames["H/y/" * string($index)] = Hyn
    frames["M/x/" * string($index)] = Mxn
    frames["M/y/" * string($index)] = Myn
    frames["phi/" * string($index)] = phin
    frames["angles/" * string($index)] = angles
    frames["F/x/" * string($index)] = fxn
    frames["F/y/" * string($index)] = fyn
    frames["t/" * string($index)] = float($index * dt)
  end
end

"""
A função `transient(...)` realiza a simulação do sistema do tempo zero até o tempo adimensional `t`.
`n` é o tamanho da malha
`dt` é o passo de tempo
`Re` é o número de Reynolds
`Cpm` é o ...
`$\alpha$` é um parâmetro parte da definição de $\mathbf{M}_0 = M_S L(\alpha|\mathbf{H}|)\hat{e}_H$, sendo $L$ a função de Langevin e $\hat{e}_H$ um vetor unitário na direção do campo magnético aplicado
`a` e `b` são a posição `x` e `y` to centro do campo magnético aplicado
`save` é uma variável booleana que indica que um arquivo com os dados de simulação deve ser salvo. Caso falso, só um arquivo com dados da evolução da magnetização é salvo, mas ele não contém as matrizes da evolução temporal completa para posterior análise
`c1` é definido por $\frac{L}{\tau U}$, sendo $\tau$ o tempo de relaxação rotacional de movimento browniano. Da adimensionalização, tem-se que L é o tamanho característico e U a velocidade característica.
`fps` quer dizer frames per second. O dt pode ser bem pequeno e pode-se não desejar salvar todos os frames, daí se escolhe uma quantidade específica por segundo adimensional para se salvar.
`filename` é o nome do arquivo hdf5 no qual serão salvos os dados de simulação
`convec` é usado para indicar se o termo convectivo deve ser considerado. Deve ser igual a 1.0 ou 0.0.
"""
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

  "Instante de início da simulação"
  println(Libc.strftime(time()), "\n")

	steps = round(Int,t/dt)

  "Índice que contém o ponto central da malha não escalonada que será impresso na tela"
	c = round(Int,n/2);

  NS = createNSObject(n, float(Re))

	for i in -1:n-2
		NS.uB[i+2] = 0*(sinpi(i*dx))^2 #inicialmente, repouso
	end

    #Condições de contorno, basta editar em v_old
	NS.v_old.x[:,n] = 2*NS.uB

    #Força
	#fx = zeros(n,n)
	#fy = zeros(n,n)

	un = zeros(n-2,n-2)+1e-15 #malha não escalonada
	vn = zeros(n-2,n-2)+1e-15
	pn = zeros(n-2,n-2)+1e-15

	numberFrames = round(Int,fps*t)
	timeToSave = round(Int,steps/numberFrames)

  #Variáveis para a parte magnética
  phi = zeros(n,n)
  Mx = zeros(n,n)
  Mx_old = zeros(n,n)
  My = zeros(n,n)
  My_old = zeros(n,n);
  left = zeros(n)
  right = zeros(n)
  upper = zeros(n)
  lower = zeros(n)
  for i in -1:n-2
    x = (i+0.5)*dx
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
   #In Julia, these variables will be in the scope out of this if clause
   file = h5open(filename, "w")
   frames = g_create(file, "frames")
   attrs(frames)["Description"] = "Frames with time series of several variables. Naming follows convention of 'variable name/direction x or y if applicable/frameNumber'. Notice direction does not apply to pressure, angles and potential field."
   frameNumber = 0
  end

  if(save)
      @saveFrames(0, dt)
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
        @saveFrames(i, dt)
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
