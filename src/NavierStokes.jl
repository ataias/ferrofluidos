module NavierStokes

#export navier_stokes_step1!, navier_stokes_step2!, navier_stokes_step3!
export solve_navier_stokes!
export staggered2not!, simpson

using Poisson
using NavierTypes

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
				u[i,j]  =  (mu*(u_s-4*u_old[i,j])/dx2+fx[i,j])*(dt/rho)
				u[i,j] += -dtx*u_old[i,j]*(u_old[i+1,j]-u_old[i-1,j])
				u[i,j] += -dtx*v_t*(u_old[i,j+1]-u_old[i,j-1])
				u[i,j] +=  u_old[i,j]
			end #if
			#Now, let's compute v in i,j
			if j!=2 && j!=n #in the case j==2, it is a boundary point
				v[i,j]  =  (mu*(v_s-4*v_old[i,j])/dx2+fy[i,j])*(dt/rho)
				v[i,j] += -dtx*u_t*(v_old[i+1,j]-v_old[i-1,j])
				v[i,j] += -dtx*v_old[i,j]*(v_old[i,j+1]-v_old[i,j-1])
				v[i,j] += v_old[i,j]
			end #if
		end #for j
	end #for i

end #navier_stokes_step1

function navier_stokes_step2!(n, A, dt, mu, rho, p, #pressão ter sido inicializada
						 	  u, v, #u e v estrela, obtido do passo 1
						 	  u_old, v_old, #u antes de u estrela
						 	  fx, fy)
	DIVij = 0.0
	dx = 1/(n-2) 
	dx2 = dx*dx
	dtx = dt/(2*dx)
	rhodtdx = rho/dt/dx
    
    #nabla^2 p = f
    f = zeros(n,n)
	for i in 2:n-1
		for j in 2:n-1
            f[i,j] = rhodtdx*((u[i+1,j]-u[i,j])+(v[i,j+1]-v[i,j]))
		end
	end

    #Condições de fronteira de fluxo
    left = zeros(n)
    right = zeros(n)
    upper = zeros(n)
    lower = zeros(n)

	#Processar fronteira esquerda
	i = 2
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i+1,j]+4*u_old[i+2,j]-u_old[i+3,j]
		left[j] = (mu/dx)*sum_aux + rho*dx*fx[i,j]
	end

	#Processar fronteira direita
	i = n #n é o tamanho da malha escalonada
	for j in 2:n-1
		sum_aux = 2*u_old[i,j]-5*u_old[i-1,j]+4*u_old[i-2,j]-u_old[i-3,j]
		right[j] = (mu/dx)*sum_aux + rho*dx*fx[i,j]
	end
	
	#Processar fronteira inferior
	j = 2
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j+1]+4*v_old[i,j+2]-v_old[i,j+3]
		lower[i] = (mu/dx)*sum_aux + rho*dx*fy[i,j]
	end

	#Processar fronteira superior
	j = n #n é o tamanho da malha escalonada
	for i in 2:n-1
		sum_aux = 2*v_old[i,j]-5*v_old[i,j-1]+4*v_old[i,j-2]-v_old[i,j-3]
		upper[i] = (mu/dx)*sum_aux + rho*dx*fy[i,j]
	end

    poissonNeumannSparseSolver(n, p, f, A, left, right, upper, lower)
    
    
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
			px = (p[i,j]-p[i-1,j])/(dx)
			py = (p[i,j]-p[i,j-1])/(dx)
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

#function solve_navier_stokes!(n, A, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)
function solve_navier_stokes!(NS::NSEquation)
        
	# ------------------------ Passo 1 ------------------------
	navier_stokes_step1!(NS.params.n, NS.params.dt, NS.params.mu, NS.params.rho, 
                         NS.v.x, NS.v.y, NS.v_old.x, NS.v_old.y, NS.f.x, NS.f.y, NS.uB)

	# ------------------------ Passo 2 ------------------------
    navier_stokes_step2!(NS.params.n, NS.msystem.chol, NS.params.dt, NS.params.mu, 
                         NS.params.rho, NS.p, NS.v.x, NS.v.y, NS.v_old.x, NS.v_old.y, 
                         NS.f.x, NS.f.y)

	# ------------------------ Passo 3 ------------------------
	navier_stokes_step3!(NS.params.n, NS.params.dt, NS.params.mu, NS.params.rho, 
                         NS.p, NS.v.x, NS.v.y, NS.v_old.x, NS.v_old.y, 
                         NS.f.x, NS.f.y, NS.uB)

end

#staggered2notS
#salva pontos internos da malha escalonada
#interpolação linear é feita para obter pontos na quina esquerda-inferior de cada bloco
#u é a entrada, malha escalonada, dimensão n
#un é a saída, dimensão n-2, escalonada de dimensão menor
#n é a dimensão da malha escalonada
function staggered2not!(u, v, p, un, vn, pn, n) 
	for i in 2:n-1 
		for j in 2:n-1 
			un[i-1,j-1] = (u[i,j] + u[i,j-1])/2
			vn[i-1,j-1] = (v[i,j] + v[i-1,j])/2
			pn[i-1,j-1] = p[i,j] #continua no centro...
		end 
	end
end


#f é um vetor com o valor a ser integrado
#este vetor tem n pontos e a distância entre eles é de 1/n
function simpson(f, n)
	#Only works for even n
    #Approximates the definite integral of f from a to b by
    #the composite Simpson's rule, using n subintervals
    h = 1 / n
    s = f[1] + f[n]
 
    for i in 2:2:n-2
        s += 4 * f[i]
    end

    for i in 3:2:n-1
        s += 2 * f[i]
    end
 
    return s * h / 3
end

end