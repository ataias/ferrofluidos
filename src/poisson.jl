module Poisson
#This is used to solve the system Ax = b in step 2
export poissonNeumannSparseSolver, getANeumannSparse, testPoissonNeumannSparse
export solveDirichletPoissonExplicit!, testPoissonDirichletExplicit
### Poisson Neumann sparse begin ----------------------------------------
function getANeumannSparse(n)
    spn = (n-2)*(n-2) #sparse n
    #Definir matriz A
    A = 4*speye(spn, spn)
    
    k=0
    for i in 2:n-1
        for j in 2:n-1
            k = j+(i-2)*(n-2) - 1
            if i-1==1 || i+1==n
                A[k,k] -= 1
            end
            
            if j-1==1 || j+1==n
                A[k,k] -= 1
            end
        end
    end
    
    for i in 1:spn-1
		if i % (n-2) != 0 
            A[i,i+1]=-1
            A[i+1,i]=-1
        end
    end
    
    for i in 1:(spn-(n-2))
		A[i,i+n-2]=-1
        A[i+n-2,i]=-1
    end
    
    
    #Escolhe-se um ponto de u e diz-se que o valor é conhecido
    #Por conveniência, escolho 0, escolho o ponto p_12 arbitrariamente
    #Nem todo ponto funciona, é bom verificar se a matriz é positiva definida
    # com o ponto escolhido!
    A[1,1] = 3 #pelo fato de eu dizer que p_12 = 0
    
    return A
end

#A can be already factorized...
function poissonNeumannSparseSolver(n, p, f, A, left, right, upper, lower)
    spn = (n-2)*(n-2) #sparse n
    dx = 1/(n-2)
    dx2 = dx*dx
    
    #Construct b vector
    b = zeros(spn) #remember, we need to solve Ax = b
    k=0
    for i in 2:n-1
        for j in 2:n-1
            k = j+(i-2)*(n-2) - 1
            b[k] = dx2*f[i,j]
            if i-1==1
                b[k] += dx*left[j]
            elseif i+1==n
                b[k] -= dx*right[j]
            end
            
            if j-1==1 
                b[k] += dx*lower[i]
            elseif j+1==n
                b[k] -= dx*upper[i]
            end
        end
    end
            
    #p_12 was considered to be 0, then:
    b[1] -= dx*left[2]
    
    #A was modified in order to be positive definite,
    #before, it was negative definite
    x = A\b;
            
    #arrange result in matrix
    for i in 2:n-1
        for j in 2:n-1
            k = j+(i-2)*(n-2) - 1
            p[i,j] = -x[k]
        end
    end
            
    #Processar fronteira esquerda
	i = 2
	for j in 2:n-1
		p[i-1,j] = p[i,j] - dx*left[j]
	end

	#Processar fronteira direita
	i = n #n é o tamanho da malha escalonada
	for j in 2:n-1
		p[i,j] = p[i-1,j] + dx*right[j]
	end

	#Processar fronteira inferior
	j = 2
	for i in 2:n-1
		p[i,j-1] = p[i,j] - dx*lower[i]
	end

	#Processar fronteira superior
	j = n #n é o tamanho da malha escalonada
	for i in 2:n-1
		p[i,j] = p[i,j-1] + dx*upper[i]
	end
       
end #poissonNeumannSparseSolver

function testPoissonNeumannSparse(n = 102)
    dx = 1/(n-2);
    p_solucao = zeros(n-2,n-2);
    for i in 1:n-2
      for j in 1:n-2
          #pressão avaliada no centro da célula
          x = (i-1)*dx + (0.5)*dx;
          y = (j-1)*dx + (0.5)*dx;
          p_solucao[i,j] = (cosh(2*pi*x)/(2*pi*sinh(2*pi)))*cos(2*pi*y);
      end
    end

    p = zeros(n,n); #valor da pressão no meio da célula
    f = zeros(n,n); #valor da f também no meio da célula
    left = zeros(1,n);
    right = zeros(1,n);
    for i in -1:n-2
      x = (i+0.5)*dx;
      right[i+2] = cos(2*pi*x);
    end
    upper = zeros(1,n);
    lower = zeros(1,n);
    
    A =getANeumannSparse(n)

    #Solve for p with Neumann conditions here
    poissonNeumannSparseSolver(n, p, f, A, left, right, upper, lower)
            
    pn = zeros(n-2, n-2);
    staggered2not!(p,pn,n)
    pn[:,:] = pn - mean(pn)
    
    println("Error = ", maximum(abs(p_solucao-pn)))
end

### Poisson Neumann sparse end ----------------------------------------

### Poisson Dirichlet explicit begin ----------------------------------------

function poissonDirichletStep!(n, p, f, left, right, upper, lower)

	dx = 1/(n-2) 
	dx2 = dx*dx
	# Processar pontos internos
	for i in 2:n-1
		for j in 2:n-1
			sum_aux = (p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1])/4
			p_new = sum_aux-dx2*f[i,j]/4
            #SOR is not working with that code!
            p[i,j] = p_new
		end
	end

	#Processar fronteira esquerda
	i = 2
	for j in 2:n-1
		p[i-1,j] = (-6*p[i,j] + p[i+1,j] + 8*left[j])/3
	end

	#Processar fronteira direita
	i = n #n é o tamanho da malha escalonada
	for j in 2:n-1
		p[i,j] = (-6*p[i-1,j] + p[i-2,j] + 8*right[j])/3
	end
	
	#Processar fronteira inferior
	j = 2
	for i in 2:n-1
		p[i,j-1] = (-6*p[i,j] + p[i,j+1] + 8*lower[i])/3
	end

	#Processar fronteira superior
	j = n #n é o tamanho da malha escalonada
	for i in 2:n-1
		p[i,j] = (-6*p[i,j-1] + p[i,j-2] + 8*upper[i])/3
	end
end

#Falta corrigir com as novas derivadas

function testPoisson(n, p, f, left, right, upper, lower)
	m = 0.0
	dx = 1/(n-2)
	dx2 = dx*dx
	for i in 2:n-1
		for j in 2:n-1
			m_aux = 0.25*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]-dx2*f[i,j])-p[i,j]
			if abs(m_aux) > m;	m = abs(m_aux); end
		end
	end
	return m
end

function solveDirichletPoissonExplicit!(n, p, f, left, right, upper, lower)
	# ------------------------ Método Iterativo ------------------------
	error = 1.0
	threshold = 1e-10
    
	i = 0
	while error > threshold
		i = i + 1 #identificar iteração
		poissonDirichletStep!(n, p, f, left, right, upper, lower)
		error = testPoisson(n, p, f, left, right, upper, lower)
		if i == 1000000 #Evita loops que sejam muito longos
			println("Problema na convergência da pressão")
			return(1)
		end
	end
end

function testPoissonDirichletExplicit(n = 52)
      dx = 1/(n-2);
      p_solucao = zeros(n-2,n-2);
      for i in 1:n-2
          for j in 1:n-2
              #pressão avaliada no centro da célula
              x = (i-1)*dx + (0.5)*dx;
              y = (j-1)*dx + (0.5)*dx;
              p_solucao[i,j] = (sinh(pi*y)/(sinh(pi)))*sin(pi*x);
          end
      end

      p_staggered = zeros(n,n); #valor da pressão no meio da célula
      f = zeros(n,n); #valor da f também no meio da célula
      left = zeros(n);
      right = zeros(n);
      upper = zeros(n);
      lower = zeros(n);
      for i in -1:n-2
          x = (i+0.5)*dx;
          upper[i+2] = sinpi(x);
      end
      

    solveDirichletPoissonExplicit!(n, p_staggered, f, left, right, upper, lower)
    pn = zeros(n-2, n-2);
    staggered2not!(p_staggered,pn,n)
    return maximum(abs(p_solucao - pn))
end

### Poisson Dirichlet explicit end ----------------------------------------

function staggered2not!(p, pn, n)
	for i in 2:n-1
		for j in 2:n-1
			pn[i-1,j-1] = p[i,j]
		end
	end
end #staggered2not!

end #module