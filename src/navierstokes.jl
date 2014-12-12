module NavierStokes

export navier_stokes_step1!, navier_stokes_step2!, navier_stokes_step3!
export solve_navier_stokes!, testPoisson, isdXok, getDt

#navier_stokes_step1!
#Obtem u* e v*, a velocidade antes de se considerar a pressão
#Note que as fronteiras ainda não 
function navier_stokes_step1!(n, dt, mu, rho, u, v, u_old, v_old, fx, fy, uB)
	#u and v will be modified!!!
	#int i, j; double u_s, v_s, u_t, v_t
	dx = 1/(n-2)
	dx2 = dx*dx
	dtx = dt/(2*dx)


	#Boundary conditions, velocidades normais na malha escalonada
	for j in 2:n-1
		u[2,j]=0 #esquerda
		u[n,j]=0 #direita
	end

	for i in 2:n-1
		v[i,2]=0 #embaixo
		v[i,n]=0 #em cima
	end

	#Process internal points
	for i in 2:n-1
        #println("i is $(i)")
		for j in 2:n-1
			# Those two lines are just summing the stencil around i,j
			u_s  = u_old[i+1,j]+u_old[i-1,j]+u_old[i,j+1]+u_old[i,j-1]
			v_s  = v_old[i+1,j]+v_old[i-1,j]+v_old[i,j+1]+v_old[i,j-1]

			#Those lines take an average of u around i,j considering the 
			#staggered grid
			u_t  = u_old[i,j]+u_old[i+1,j]+u_old[i+1,j-1]+u_old[i,j-1]
			u_t *= 0.25
			v_t  = v_old[i,j]+v_old[i-1,j]+v_old[i-1,j+1]+v_old[i,j+1]
			v_t *= 0.25

			#Now, let's compute u in i,j
			#value is fixed in the boundaries
			if i!=2 && i!=n #in the case i==2 or n, it is a boundary point
				#Order is important here!
				u[i,j]  =  mu*(u_s-4*u_old[i,j])/dx2+fx[i,j]
				u[i,j] *= dt/rho
				u[i,j] += -dtx*u_old[i,j]*(u_old[i+1,j]-u_old[i-1,j])
				u[i,j] += -dtx*v_t*(u_old[i,j+1]-u_old[i,j-1])
				u[i,j] +=  u_old[i,j]
			end #if
			#Now, let's compute v in i,j
			if j!=2 && j!=n #in the case j==2 or n, it is a boundary point
				v[i,j]  =  mu*(v_s-4*v_old[i,j])/dx2+fy[i,j]
				v[i,j] *= dt/rho
				v[i,j] += -dtx*u_t*(v_old[i+1,j]-v_old[i-1,j])
				v[i,j] += -dtx*v_old[i,j]*(v_old[i,j+1]-v_old[i,j-1])
				v[i,j] += u_old[i,j]
			end #if
		end #for j
	end #for i
    println("u=", u)
    println("v=", v)
end #navier_stokes_step1

function navier_stokes_step2!(n, dt, mu, rho, p, #pressão pode ser matriz zeros(6,6)
						 	  u, v, #u e v estrela, obtido do passo 1
						 	  u_old, v_old, #u antes de u estrela
						 	  fx, fy)
	DIVij = 0.0
	dx = 1.0/(n-2.0) #n-1 or n-2?
	dx2 = dx*dx
	dtx = dt/(2.*dx)

	# Processar pontos internos
	for i in 2:n-1
		for j in 2:n-1
			DIVij = (rho/dt)*(u[i+1,j]-u[i,j]+v[i,j+1]-v[i,j])/dx
			sum_aux = p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]
			p[i,j] = 0.25*sum_aux-0.25*dx2*DIVij
		end
	end

	#Processar fronteira i=2
	i = 2
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i+1,j]+4*u_old[i+2,j]-u_old[i+3,j]
		f_aux = fx[i,j]+fx[i+1,j]
		p[i-1,j] = p[i,j]-(mu/dx)*sum_aux-(rho*dx/2)*f_aux
	end

	#Processar fronteira i=n-1
	i = n-1 #n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i-1,j]+4*u_old[i-2,j]-u_old[i-3,j]
		f_aux = fx[i,j]+fx[i-1,j]
		p[i+1,j] = p[i,j]+(mu/dx)*sum_aux+(rho*dx/2)*f_aux
	end
	
	#Processar fronteira j=2
	j = 2
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j+1]+4*v_old[i,j+2]-v_old[i,j+3]
		f_aux = fy[i,j]+fy[i,j+1]
		p[i,j-1] = p[i,j]-(mu/dx)*sum_aux-(rho*dx/2)*f_aux
	end

	#Processar fronteira j=n-1
	j = n-1 #n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j-1]+4*v_old[i,j-2]-v_old[i,j-3]
		f_aux = fy[i,j]+fy[i,j-1]
		p[i,j+1] = p[i,j]+(mu/dx)*sum_aux+(rho*dx/2)*f_aux
	end
end

function navier_stokes_step3!(n, dt, mu, rho, p, #pressão já calculada
						 	  u, v, #u* e v*, serão atualizados para u(n+1)
						 	  u_old, v_old, #u(n)
						 	  fx, fy, uB)
	dx = 1/(n-2)
	px = 0.0
	py = 0.0
	#Pontos internos
	for i in 2:n-1
		for j in 2:n-1
			px = (p[i,j]-p[i-1,j])/dx
			py = (p[i,j]-p[i,j-1])/dx
			u[i,j] = u[i,j] - (dt/rho)*px
			v[i,j] = v[i,j] - (dt/rho)*py
		end
	end

	#Fronteira
	for j in 2:n-1 v[1,j] = -v[2,j] end #esquerda
	for i in 2:n-1 u[i,1] = -u[i,2] end #embaixo
	for j in 2:n-1 v[n,j] = -v[n-1,j] end #direita
	for i in 2:n-1 u[i,n] = 2*uB[i]-u[i,n-1] end # em cima
end

function solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)
	# ------------------------ Passo 1 ------------------------
	navier_stokes_step1!(n, dt, mu, rho, u, v, u_old, v_old, fx, fy, uB)

	# ------------------------ Passo 2 ------------------------
	#A parte 2 tem de ter iterações até convergir
	p_old = zeros(n,n)
	error = 1.0
	threshold = 1e-8
	i = 0
	while error > threshold
		navier_stokes_step2!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy)
		p[:,:] = p - minimum(p)
		p[1,1] = 0; p[1,n]=0; p[n,1]=0; p[n,n]=0
		error = maxabs(p-p_old)

		p_old[:,:] = p
		i = i + 1
		if i == 12000 #limitado a 12 mil iterações... pode mudar este número
			return "Fim"
		end
	end
	println("There were $(i) iterations to solve step2. Estimated Error=$(error)")
	println("Actual error: $(testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy))")

	# ------------------------ Passo 3 ------------------------
	navier_stokes_step3!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)

end

#testPoisson
#
function testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy)
	m = zeros(n,n)
	dx = 1/(n-2)
	dx2 = dx^2
	for i in 2:n-1
		for j in 2:n-1
			DIVij = (rho/dt)*(u[i+1,j]-u[i,j]+v[i,j+1]-v[i,j])/dx #este é o u e v estrela
			m[i,j] = 0.25*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]-dx2*DIVij)-p[i,j]
		end
	end

	#Processar fronteira i=2
	i = 2
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i+1,j]+4*u_old[i+2,j]-u_old[i+3,j]
		f_aux = fx[i,j]+fx[i+1,j]
		m[i-1,j] = p[i,j]-(mu/dx)*sum_aux-(rho*dx/2)*f_aux - p[i-1,j]
	end

	#Processar fronteira i=n-1
	i = n-1 #n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i-1,j]+4*u_old[i-2,j]-u_old[i-3,j]
		f_aux = fx[i,j]+fx[i-1,j]
		m[i+1,j] = p[i,j]+(mu/dx)*sum_aux+(rho*dx/2)*f_aux - p[i+1,j]
	end
	
	#Processar fronteira j=2
	j = 2
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j+1]+4*v_old[i,j+2]-v_old[i,j+3]
		f_aux = fy[i,j]+fy[i,j+1]
		m[i,j-1] = p[i,j]-(mu/dx)*sum_aux-(rho*dx/2)*f_aux - p[i,j-1]
	end

	#Processar fronteira j=n-1
	j = n-1 #n é o tamanho da malha escalonada
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j-1]+4*v_old[i,j-2]-v_old[i,j-3]
		f_aux = fy[i,j]+fy[i,j-1]
		m[i,j+1] = p[i,j]+(mu/dx)*sum_aux+(rho*dx/2)*f_aux - p[i,j+1]
	end
	println(m)
	return maxabs(m)
end

#isdXok
#Esta função analisa se o número de pontos escolhido 
#satisfaz a condição de camada limite hidrodinâmica
function isdXok(Re, n)
	dx = 1./(n-2); #because n = n' + 1
	if dx < (1/Re) 
		return true
	else 
		return false
	end
end

function getDt(n, Re, divFactor=5)
	dx = 1/(n-2) #n é o tamanho da malha escalonada, por isso está n-2 ao invés de n-1
	dt1 = 0.25*dx*dx/Re
	dt2 = dx
	dt = 0.0
	if dt1 < dt2
		dt = dt1/divFactor
	else 
		dt = dt2/divFactor
	end
	println("dt = ", dt);
	return dt;
end

#staggered2notS
#descarta dimensão extra da malha escalonada
#u é a entrada, malha escalonada, dimensão n
#un é a saída, dimensão n-1, escalonada de dimensão menor
#n é a dimensão da malha escalonada
function staggered2notS!(u, un, n) 
	for i in 2:n 
		for j in 2:n un[i-1,j-1] = u[i,j] end 
	end
end
end