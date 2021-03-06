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
    return i <= 100 ? i/100.0 : 1.0
end

macro saveFrames(i, time)
  return quote
    frames["velocity/x/" * string($i)] = un
    frames["velocity/y/" * string($i)] = vn
    frames["pressure/" * string($i)] = pn
    frames["H/x/" * string($i)] = Hxn
    frames["H/y/" * string($i)] = Hyn
    frames["M/x/" * string($i)] = Mxn
    frames["M/y/" * string($i)] = Myn
    frames["phi/" * string($i)] = phin
    frames["angles/" * string($i)] = angles
    frames["F/x/" * string($i)] = fxn
    frames["F/y/" * string($i)] = fyn
    frames["t/" * string($i)] = float($time)

    #Trying to change an existing value gives an error message
    if exists(frames, "lastFrame")
      o_delete(frames, "lastFrame")
    end
    frames["lastFrame"] = $i
  end
end

"""
A função `transient(...)` realiza a simulação do sistema do tempo zero até o tempo adimensional `t`.
`n` é o tamanho da malha
`dt` é o passo de tempo
`Re` é o número de Reynolds
`Cpm` é o ...
`\alpha` é um parâmetro parte da definição de \mathbf{M}_0 = M_S L(\alpha|\mathbf{H}|)\hat{e}_H, sendo L a função de Langevin e \hat{e}_H um vetor unitário na direção do campo magnético aplicado
`a` e `b` são a posição `x` e `y` to centro do campo magnético aplicado
`save` é uma variável booleana que indica que um arquivo com os dados de simulação deve ser salvo. Caso falso, só um arquivo com dados da evolução da magnetização é salvo, mas ele não contém as matrizes da evolução temporal completa para posterior análise
`c1` é definido por L/(\tau U), sendo \tau o tempo de relaxação rotacional de movimento browniano. Da adimensionalização, tem-se que L é o tamanho característico e U a velocidade característica.
`fps` quer dizer frames per second. O dt pode ser bem pequeno e pode-se não desejar salvar todos os frames, daí se escolhe uma quantidade específica por segundo adimensional para se salvar.
`filename` é o nome do arquivo hdf5 no qual serão salvos os dados de simulação
`convec` é usado para indicar se o termo convectivo deve ser considerado. Deve ser igual a 1.0 ou 0.0.
"""
function transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, filename, convec=0.0; should_print=true)

  dx = 1/(n-2)

  if should_print
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
    println(" convec\t= ", convec)
  end
  date = Libc.strftime(time())

  if(save)
    file = h5open(filename, "w")
    info = g_create(file, "info")
    attrs(info)["Description"] = "Info about simulation, like date and time of simulation."
    info["date"] = date
    info["n"] = n - 2
    info["dx"] = dx
    info["t"] = t
    info["dt"] = dt
    info["Re"] = Re
    info["Cpm"] = Cpm
    info["alpha"] = alpha
    info["a"] = a
    info["b"] = b
    info["c1"] = c1
    info["Pe"] = 1/c1
    info["fps"] = fps
    info["completed"] = 0
  end

  if should_print
    "Instante de início da simulação"
    println(date, "\n")
  end

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

  # A small offset is added to avoid complaining from matplotlib if it there
  # are only zeros
  fact = 0
  angles = zeros(n-2, n-2)+1e-15;
  fxn = zeros(n-2, n-2)+1e-15
  fyn = zeros(n-2, n-2)+1e-15

  if(save)
   #In Julia, these variables will be in the scope out of this if clause
   frames = g_create(file, "frames")
   attrs(frames)["Description"] = "Frames with time series of several variables. Naming follows convention of 'variable name/direction x or y if applicable/frameNumber'. Notice direction does not apply to pressure, angles and potential field."
  end
  frameNumber = 0

  if(save)
    #Salvar valores das matrizes em A
    @saveFrames(0, 0.0)
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

  # Auxiliar variable to compute force
  tau = zeros(n-2)

  # TODO: Create variable to store vorticity, it can be in memory
  # 24/05/2016
  midPoint = Dict(:t     => zeros(numberFrames),
                  :vortc => zeros(numberFrames),
                  :u     => zeros(numberFrames),
                  :v     => zeros(numberFrames))

	for i in 1:steps
    fact = factor(i)
    for j in -1:n-2
        NS.uB[j+2] = fact*(sinpi(j*dx))^2
    end

    #Obter força
    #O primeiro passo é obter M em cada passo de tempo
    #faz sentido não usar o fator fact para getM, pois a evolução inicia
    #do valor anterior de M, que é 0 no tempo 0
    if Cpm > 1e-8
      v∇M!(n, NS.v.x, NS.v.y, Mx_old, My_old, v∇Mx, v∇My)
      getM!(n, c1, dt, Mx, My, Mx_old, My_old, Mx0, My0, convec*v∇Mx, convec*v∇My)
      getPhi!(n, phi, Mx, My, fHx, fHy, A)
      getH!(n, phi, Hx, Hy)
      getForce!(n, Cpm, Hx, Hy, Mx, My, NS.f.x, NS.f.y)
    end
    # Resolver hidrodinâmica
    solve_navier_stokes!(NS)

		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(NS.v.x, NS.v.y, NS.p, un,  vn,  pn,   n)
      staggered2not!(Hx,     Hy,     phi,  Hxn, Hyn, phin, n)
      staggered2not!(Mx,     My,           Mxn, Myn,       n)
      staggered2not!(NS.f.x, NS.f.y,       fxn, fyn,       n)
      Angle!(Hxn, Hyn, Mxn, Myn, angles, n-2)
      frameNumber += 1
      if(save)
        @saveFrames(frameNumber, i*dt)
      end #if(save)

      time = i*dt
      #Observe que vn e un estão salvas nas quinas.
      #Este cálculo de vorticidade só é válido caso n seja PAR!
			vortc =  ((vn[c+1,c]-vn[c-1,c]) - (un[c,c+1]-un[c,c-1]) )/(2*dx)
			for i in 2:n-1
				tau[i-1] = (1/Re)*((NS.v.x[i,n]-NS.v.x[i,n-1])/dx)
			end #for i in 2:n-1

			F = simpson(tau, n-2)

      if should_print
        # Presenting information
        println("t = ", time)
  			println("  u[0.5,0.5]\t= ", un[c,c])
  			println("  v[0.5,0.5]\t= ", vn[c,c])
  			println("  ω [0.5,0.5]\t= ", vortc)
        println("  Pressure values in range\t [", minimum(pn), ", ", maximum(pn), "]")
  			println("  F\t= ", F)
      end

      # Testes unitários -> retorno da função será valores do ponto intermediário
      # midPoint[:t][frameNumber] = time
      # midPoint[:u][frameNumber] = un[c,c]
      # midPoint[:v][frameNumber] = vn[c,c]
      # midPoint[:vortc][frameNumber] = vortc
		end #if (i % timeToSave == 0) || (i == steps)

    #Preparing for next time step
		NS.v.x, NS.v_old.x = NS.v_old.x, NS.v.x
		NS.v.y, NS.v_old.y = NS.v_old.y, NS.v.y

    Mx, Mx_old = Mx_old, Mx
    My, My_old = My_old, My

	end #for i in 1:steps

  if(save)
    if exists(info, "completed")
      o_delete(info, "completed")
    end
    info["completed"] = 1
    close(file)
    println("File closed")
  end

	return midPoint
end #function transient

end #module
