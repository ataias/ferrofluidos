using NavierStokes

#O mesmo dt será usado para todos
Re = 10.0;
dt = getDt(81, Re, 1.01); #escolho o menor dt entre todos	
t = 2.5;

solution = steadyState(81, dt, Re, t)
fileErro = open("errorUVvsDx2.dat", "w")

for n in 13:4:69
	dx = 1/(n-2)
	dx2 = dx^2
	ss = steadyState(n, dt, Re, t)
	write(fileErro, dx2)
	write(fileErro, abs(solution[1] - ss[1]))
	write(fileErro, abs(solution[2] - ss[2]))
end

close(fileErro)

#retorna o valor em regime permanente do ponto 0.5, 0.5
#resolve equações para um dado n e Re
#t is time, in seconds, of physical simulation
function steadyState(n, dt, Re, t)
	println("Running! n = ", n)
	rho = 1.0;
	mu = 1/Re;
	dx = 1/(n-2)
	steps = integer(t/dt)

	if(!isdXok(Re, n))
		println("There is a problem with your dx. Increase n.\n");
		return 0
	end

	p = zeros(n,n);
	u = zeros(n,n);
	uB = zeros(1,n);
	for i in -1:n-2
		uB[i+2] = (sinpi(i*dx))^2
	end
	u[:,n] = 2*uB
	u[1,:] = NaN
		
	v = zeros(n,n);
	v[:,1] = NaN

	u_old = copy(u);
	v_old = copy(v);
	fx = zeros(n,n);
	fy = zeros(n,n);

	# file = open("N" * string(n-2) * ".dat", "w")
	#falta ainda salvar o valor do ponto médio conforme
	#o tempo avança...

	un = zeros(n-2,n-2)
	vn = zeros(n-2,n-2)
	pn = zeros(n-2,n-2)

	numberFrames = integer(180*t)
	timeToSave = integer(steps/numberFrames)

	for i in 1:steps
		#println("Step i = ", i)
		solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)

		if i % timeToSave == 0
		 	println("t = ", i*dt)
		# 	staggered2not!(u, v, p, un, vn, pn, n)
		# 	write(file, un) #note that it writes to memory column-wise
		# 	write(file, vn)
		# 	write(file, pn)
		end

		u, u_old = u_old, u
		u[:,:] = u_old[:,:]
	end

	# close(file)

	#steady state value in (0.5, 0.5)
	staggered2not!(u, v, p, un, vn, pn, n)
	ij = integer(n/2 - 1)
	return [un[ij, ij], vn[ij, ij]] #returns steady state value in the middle
end