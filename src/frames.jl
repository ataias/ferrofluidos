using NavierStokes

#modo de usar
#n  - ARGS[1] é o tamanho da matriz escalonada
#t  - ARGS[2] é o tempo de simulação, em segundos
#Re - ARGS[3] é o número de Reynolds
#dt - ARGS[4] dividir o passo de tempo, deve ser maior do que 1

#Exemplo:
# 		julia frames.jl 43 2.5 10.0 1.25

n = int(ARGS[1]);
t = float(ARGS[2]);
Re = float(ARGS[3]);
dt = getDt(n, Re, float(ARGS[4]));

println("Dados sobre simulação:\n n\t= ", n, "\n t\t= ", t, "\n Re\t= ", Re, "\n dt\t= ", dt, "\n");

#steadyState
#retorna o valor em regime permanente do ponto 0.5, 0.5
#resolve equações para um dado n e Re
#t is time, in seconds, of physical simulation
function steadyState(n, dt, Re, t)
	rho = 1.0
	mu = 1/Re
	dx = 1/(n-2)
	steps = integer(t/dt)
	c = integer(n/2); #center, non-staggered grid
#	ce = integer(n/2) #center, staggered grid

	if(!isdXok(Re, n))
		println("There is a problem with your dx. Increase n.\n")
		return 1
	end

	#Declarando vetores e matrizes que serão utilizados
	p = zeros(n,n)

	uB = zeros(1,n)
	for i in -1:n-2
		uB[i+2] = (sinpi(i*dx))^2
	end

	u = zeros(n,n) #malha escalonada
	u[:,n] = 2*uB
	u[1,:] = NaN
		
	v = zeros(n,n); #malha escalonada
	v[:,1] = NaN

	u_old = copy(u)
	v_old = copy(v)

	fx = zeros(n,n)
	fy = zeros(n,n)

	file = open("N" * string(n-2) * ".dat", "w")

	un = zeros(n-2,n-2) #malha não escalonada
	vn = zeros(n-2,n-2)
	pn = zeros(n-2,n-2)

	numberFrames = integer(60*t)
	timeToSave = integer(steps/numberFrames)

	for i in 1:steps
		solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)

		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(u, v, p, un, vn, pn, n)
			write(file, un); write(file, vn);
			write(file, pn)
			println("t = ", i*dt)
			println("  u[0.5,0.5]\t= ", un[c,c])
			println("  v[0.5,0.5]\t= ", vn[c,c])
			vortc =  ((vn[c+1,c]-vn[c-1,c]) - (un[c,c+1]-un[c,c-1]))/(2*dx)
			println("  ω [0.5,0.5]\t= ", vortc)
#			vortce  = ((v[ce+1,ce+1]+v[ce+1,ce]) - (v[ce-1,ce+1]+v[ce-1,ce]))/(4*dx)
#            vortce -= ((u[ce+1,ce+1]+u[ce,ce+1]) - (u[ce+1,ce-1]+u[ce,ce-1]))/(4*dx)
#			println("  ω [0.5,0.5]\t= ", vortce)
#			F = 0.0
#			for i in 2:n-1
#				F = F + (u[i,n]-u[i,n-1])/dx
#			end
#			F = F/(n-2)
#			println("  F\t= ", F)
		end

		u, u_old = u_old, u
		v, v_old = v_old, v
		u[:,:] = u_old[:,:]
		v[:,:] = v_old[:,:]
	end

	close(file)
	return 0
end

@time steadyState(n, dt, Re, t)
