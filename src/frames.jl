using NavierStokes
using NavierTypes

#modo de usar
#n   - ARGS[1] é o tamanho da matriz escalonada
#t   - ARGS[2] é o tempo de simulação, em segundos
#Re  - ARGS[3] é o número de Reynolds
#dt  - ARGS[4] dividir o passo de tempo, deve ser maior do que 1
#save - ARGS[5] se for 1, salva um arquivo com as matrizes evoluindo no tempo
#Exemplo:
# 		julia frames.jl 42 2.5 10.0 1.25 0
# Se quiser salvar num arquivo a saída do terminal:
#       julia frames.jl 42 2.5 10.0 1.25 0 >> out.txt

n = int(ARGS[1]);
t = float(ARGS[2]);
Re = float(ARGS[3]);
divFactor = float(ARGS[4]);
dt = getDt(n, Re, float(ARGS[4]));
save = bool(int(ARGS[5]))

println("Dados sobre simulação:\n n\t= ", n, "\n dx\t= ", 1/(n-2), "\n t\t= ", t, "\n Re\t= ", Re, "\n dt\t= ", dt, "\n ", strftime(time()), "\n");

#steadyState
#retorna o valor em regime permanente do ponto 0.5, 0.5
#resolve equações para um dado n e Re
#t is time, in seconds, of physical simulation
function steadyState(n, dt, Re, t, save)

	dx = 1/(n-2)
	steps = integer(t/dt)
	c = integer(n/2); #center, non-staggered grid
#	ce = integer(n/2) #center, staggered grid

    NS = createNSObject(n, Re, divFactor)

	for i in -1:n-2
		NS.uB[i+2] = (sinpi(i*dx))^2
	end

    #Condições de contorno, basta editar em v_old
	NS.v_old.x[:,n] = 2*NS.uB

    #Força
	#fx = zeros(n,n)
	#fy = zeros(n,n)

    if(save)
	   file = open("N" * string(n-2) * ".dat", "w")
    end

	un = zeros(n-2,n-2) #malha não escalonada
	vn = zeros(n-2,n-2)
	pn = zeros(n-2,n-2)

	numberFrames = integer(60*t)
	timeToSave = integer(steps/numberFrames)
    
	for i in 1:steps
		solve_navier_stokes!(NS)

		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(NS.v.x, NS.v.y, NS.p, un, vn, pn, n)
            if(save)
                write(file, un); write(file, vn);
                write(file, pn)
            end
			println("t = ", i*dt)
			println("  u[0.5,0.5]\t= ", un[c,c])
			println("  v[0.5,0.5]\t= ", vn[c,c])
			vortc =  ((vn[c+1,c]-vn[c-1,c]) - (un[c,c+1]-un[c,c-1]))/(2*dx)
			println("  ω [0.5,0.5]\t= ", vortc)
			tau = zeros(n-2)
			for i in 2:n-1
				tau[i-1] = (1/Re)*((NS.v.x[i,n]-NS.v.x[i,n-1])/dx)
			end
			F = simpson(tau, n-2)
			println("  F\t= ", F)
		end

        #Preparing for next time step
		NS.v.x, NS.v_old.x = NS.v_old.x, NS.v.x
		NS.v.y, NS.v_old.y = NS.v_old.y, NS.v.y
        
	end
    
    if(save)
        close(file)
    end 
	return 0
end

@time steadyState(n, dt, Re, t, save)
