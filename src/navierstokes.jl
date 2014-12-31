module NavierStokes

export navier_stokes_step1!, navier_stokes_step2!, navier_stokes_step3!
export solve_navier_stokes!, testPoisson, isdXok, getDt
export staggered2not!, simpson

#navier_stokes_step1!
#Obtem u* e v*, a velocidade antes de se considerar a pressão
#Note que as fronteiras ainda não 
function navier_stokes_step1!(n, dt, mu, rho, u, v, u_old, v_old, fx, fy, uB)
	#u and v will be modified!!!

	dx = 1/(n-2)
	dx2 = dx*dx
	dtx = dt/(2*dx)

	#Process internal points
	for i in 2:n-1
		for j in 2:n-1

			# Those two lines are just summing the stencil around i,j
			u_s  = u_old[i+1,j]+u_old[i-1,j]+u_old[i,j+1]+u_old[i,j-1]
			v_s  = v_old[i+1,j]+v_old[i-1,j]+v_old[i,j+1]+v_old[i,j-1]

			#Those lines take an average of u around i,j considering the 
			#staggered grid
			u_t  = (u_old[i,j]+u_old[i+1,j]+u_old[i+1,j-1]+u_old[i,j-1])/4
			v_t  = (v_old[i,j]+v_old[i-1,j]+v_old[i-1,j+1]+v_old[i,j+1])/4

			#Now, let's compute u in i,j
			#value is fixed in the boundaries
			if i!=2 && i!=n #in the case i==2 or n, it is a boundary point
				#Order is important here!
				u[i,j]  =  mu*(u_s-4*u_old[i,j])/dx2+fx[i,j]
				u[i,j] *= (dt/rho)
				u[i,j] += -dtx*u_old[i,j]*(u_old[i+1,j]-u_old[i-1,j])
				u[i,j] += -dtx*v_t*(u_old[i,j+1]-u_old[i,j-1])
				u[i,j] +=  u_old[i,j]
			end #if
			#Now, let's compute v in i,j
			if j!=2 && j!=n #in the case j==2, it is a boundary point
				v[i,j]  =  mu*(v_s-4*v_old[i,j])/dx2+fy[i,j]
				v[i,j] *= (dt/rho)
				v[i,j] += -dtx*u_t*(v_old[i+1,j]-v_old[i-1,j])
				v[i,j] += -dtx*v_old[i,j]*(v_old[i,j+1]-v_old[i,j-1])
				v[i,j] += v_old[i,j]
			end #if
		end #for j
	end #for i

    #Boundary conditions, velocidades normais na malha escalonada
	for j in 2:n-1
		u[2,j] = 0 #esquerda
		u[n,j] = 0 #direita
	end

	for i in 2:n-1
		v[i,2] = 0 #embaixo
		v[i,n] = 0 #em cima
	end

end #navier_stokes_step1

function navier_stokes_step2!(n, dt, mu, rho, p, #pressão pode ser matriz zeros(6,6)
						 	  u, v, #u e v estrela, obtido do passo 1
						 	  u_old, v_old, #u antes de u estrela
						 	  fx, fy)
	DIVij = 0.0
	dx = 1/(n-2) 
	dx2 = dx*dx
	dtx = dt/(2*dx)
	rhodt = rho/dt
	r = 2/(1+pi/(n-2)) #SOR constant

	# Processar pontos internos
	for i in 2:n-1
		for j in 2:n-1
			DIVij = rhodt*((u[i+1,j]-u[i,j])+(v[i,j+1]-v[i,j]))/(4*dx)
			sum_aux = (p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1])/4
			p_new = sum_aux-dx2*DIVij
			p[i,j] = (1-r)*p[i,j]+r*p_new
		end
	end

	#Processar fronteira esquerda
	i = 2
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i+1,j]+4*u_old[i+2,j]-u_old[i+3,j]
		p[i-1,j] = p[i,j] - (mu/dx)*sum_aux - rho*dx*fx[i,j]
	end

	#Processar fronteira direita
	i = n #n é o tamanho da malha escalonada
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i-1,j]+4*u_old[i-2,j]-u_old[i-3,j]
		p[i,j] = p[i-1,j] + (mu/dx)*sum_aux + rho*dx*fx[i,j]
	end
	
	#Processar fronteira inferior
	j = 2
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j+1]+4*v_old[i,j+2]-v_old[i,j+3]
		p[i,j-1] = p[i,j] - (mu/dx)*sum_aux - rho*dx*fy[i,j]
	end

	#Processar fronteira superior
	j = n #n é o tamanho da malha escalonada
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j-1]+4*v_old[i,j-2]-v_old[i,j-3]
		p[i,j] = p[i,j-1] + (mu/dx)*sum_aux + rho*dx*fy[i,j]
	end
end

function navier_stokes_step3!(n, dt, mu, rho, p, #pressão já calculada
						 	  u, v, #u* e v*, serão atualizados para u(n+1)
						 	  u_old, v_old, #u(n)
						 	  fx, fy, uB)
	dx = 1/(n-2)
	px = 0.0
	py = 0.0
	drho = dt/rho
	#Pontos internos
	for i in 2:n-1
		for j in 2:n-1
			px = (p[i,j]-p[i-1,j])/dx
			py = (p[i,j]-p[i,j-1])/dx
			u[i,j] = u[i,j] - drho*px
			v[i,j] = v[i,j] - drho*py
		end
	end
    
    #Boundary conditions, velocidades normais na malha escalonada
	for j in 2:n-1
		u[2,j] = 0 #esquerda
		u[n,j] = 0 #direita
	end

	for i in 2:n-1
		v[i,2] = 0 #embaixo
		v[i,n] = 0 #em cima
	end

	#Fronteira
	for j in 2:n v[1,j] = -v[2,j]          end #esquerda
	for i in 2:n u[i,1] = -u[i,2]          end #embaixo
	for j in 2:n v[n,j] = -v[n-1,j]        end #direita
	for i in 2:n u[i,n] = 2*uB[i]-u[i,n-1] end # em cima


end

function solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)
	# ------------------------ Passo 1 ------------------------
	navier_stokes_step1!(n, dt, mu, rho, u, v, u_old, v_old, fx, fy, uB)

	# ------------------------ Passo 2 ------------------------
	#A parte 2 tem de ter iterações até convergir
	error = 1.0
	threshold = 1e-14
	i = 0
	while error > threshold
		i = i + 1 #identificar iteração
		navier_stokes_step2!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy)
		error = testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy)

		if i == 80000 #Evita loops que sejam muito longos
			println("Fim")
			return 1
		end
	end

	p[:,:] = p - minimum(p)
	# p[:,:] = p - mean(p)
	# p[1,1] = 0; p[1,n]=0; p[n,1]=0; p[n,n]=0

	# if testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy) > (threshold*10)
	# 	println("Problem!")
	# 	println("There were $(i) iterations to solve step2. Estimated Error=$(error)")
	# 	println("Actual error: $(testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy))")
	# end
	# ------------------------ Passo 3 ------------------------
	navier_stokes_step3!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)

end

#testPoisson
#
function testPoisson(n, dt, mu, p, rho, u, v, u_old, v_old, fx, fy)
	m = 0.0
	dx = 1/(n-2)
	dx2 = dx*dx
	for i in 2:n-1
		for j in 2:n-1
			DIVij = (rho/dt)*(u[i+1,j]-u[i,j]+v[i,j+1]-v[i,j])/dx #este é o u e v estrela
			m_aux = 0.25*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]-dx2*DIVij)-p[i,j]
			if abs(m_aux) > m;	m = abs(m_aux); end
		end
	end

	#Processar fronteira esquerda
	i = 2
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i+1,j]+4*u_old[i+2,j]-u_old[i+3,j]
		m_aux = - p[i-1,j] + p[i,j] - (mu/dx)*sum_aux - rho*dx*fx[i,j]
		if abs(m_aux) > m;	m = abs(m_aux); end
	end

	#Processar fronteira direita
	i = n #n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i-1,j]+4*u_old[i-2,j]-u_old[i-3,j]
		m_aux  = - p[i,j] + p[i-1,j] + (mu/dx)*sum_aux + rho*dx*fx[i,j]
		if abs(m_aux) > m;	m = abs(m_aux); end
	end
	
	#Processar fronteira inferior
	j = 2
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j+1]+4*v_old[i,j+2]-v_old[i,j+3]
		m_aux = - p[i,j-1] + p[i,j] - (mu/dx)*sum_aux - rho*dx*fy[i,j]
		if abs(m_aux) > m;	m = abs(m_aux); end
	end

	#Processar fronteira superior
	j = n #n é o tamanho da malha escalonada
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j-1]+4*v_old[i,j-2]-v_old[i,j-3]
		m_aux = - p[i,j] + p[i,j-1] + (mu/dx)*sum_aux + rho*dx*fy[i,j]
		if abs(m_aux) > m;	m = abs(m_aux); end
	end

	return m
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
	dx = 1./(n-2) #n é o tamanho da malha escalonada, por isso está n-2 ao invés de n-1
	dt1 = 0.25*Re*dx*dx
	dt2 = dx
	dt = 0.0
	if dt1 < dt2
		dt = dt1/divFactor
	else 
		dt = dt2/divFactor
	end
	return dt;
end

#staggered2notS
#salva pontos internos da malha escalonada
#interpolação linear é feita para obter pontos no centro de cada bloco
#u é a entrada, malha escalonada, dimensão n
#un é a saída, dimensão n-2, escalonada de dimensão menor
#n é a dimensão da malha escalonada
function staggered2not!(u, v, p, un, vn, pn, n) 
	for i in 2:n-1 
		for j in 2:n-1 
			un[i-1,j-1] = (u[i,j] + u[i+1,j])/2
			vn[i-1,j-1] = (v[i,j] + v[i,j+1])/2
			pn[i-1,j-1] = p[i,j]
		end 
	end
end

function simpson(f, a, b, n)
	#Only works for even n
    #Approximates the definite integral of f from a to b by
    #the composite Simpson's rule, using n subintervals
    h = (b - a) / n
    s = f(a) + f(b)
 
    for i in 1:2:n
        s += 4 * f(a + i * h)
    end

    for i in 2:2:n-1
        s += 2 * f(a + i * h)
    end
 
    return s * h / 3
end

end