#Condições de Dirichlet

function poissonStep!(n, p, f, left, right, upper, lower)

	dx = 1/(n-2) 
	dx2 = dx*dx
	#r = 2/(1+pi/(n-2)) #SOR constant

	# Processar pontos internos
	for i in 2:n-1
		for j in 2:n-1
			sum_aux = (p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1])/4
			p_new = sum_aux-dx2*f[i,j]/4
#			p[i,j] = (1-r)*p[i,j]+r*p_new
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

#function testPoisson(n, p, f, left, right, upper, lower)
#	m = 0.0
#	dx = 1/(n-2)
#	dx2 = dx*dx
#	for i in 2:n-1
#		for j in 2:n-1
#			m_aux = 0.25*(p[i+1,j]+p[i-1,j]+p[i,j+1]+p[i,j-1]-dx2*f[i,j])-p[i,j]
#			if abs(m_aux) > m;	m = abs(m_aux); end
#		end
#	end
#
#	#Processar fronteira esquerda
#	i = 2
#	for j in 2:n-1
#		m_aux = p[i,j] + p[i-1,j] - 2*left[j]
#		if abs(m_aux) > m;	m = abs(m_aux); end
#	end
#
#	#Processar fronteira direita
#	i = n #n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
#	for j in 2:n-1
#		m_aux  = p[i,j] + p[i-1,j] - 2*right[j]
#		if abs(m_aux) > m;	m = abs(m_aux); end
#	end
#	
#	#Processar fronteira inferior
#	j = 2
#	for i in 2:n-1
#		m_aux = p[i,j-1] + p[i,j] - 2*lower[i]
#		if abs(m_aux) > m;	m = abs(m_aux); end
#	end
#
#	#Processar fronteira superior
#	j = n #n é o tamanho da malha escalonada
#	for i in 2:n-1
#		m_aux = p[i,j] + p[i,j-1] - 2*upper[i]
#		if abs(m_aux) > m;	m = abs(m_aux); end
#	end
#
#	return m
#end

function solvePoisson!(n, p, f, left, right, upper, lower)
	# ------------------------ Método Iterativo ------------------------
	error = 1.0
	threshold = 1e-15
    
	i = 0
	while error > threshold
		i = i + 1 #identificar iteração
		poissonStep!(n, p, f, left, right, upper, lower)
#		error = testPoisson(n, p, f, left, right, upper, lower)
##        println(error)
#		if i == 70000 #Evita loops que sejam muito longos
#			println("Problema na convergência da pressão")
#			return(1)
#		end
        if i == 500000
            error = 1e-16
        end
	end
    #p[:,:] = p - mean(p)

end

function problem1(n = 8)
      dx = 1/(n-2);
      p_solucao = zeros(n-2,n-2);
      triplets = fill([0.0, 0.0, 0.0], n-2, n-2);
      for i in 1:n-2
          for j in 1:n-2
              #pressão avaliada no centro da célula
              x = (i-1)*dx + (0.5)*dx;
              y = (j-1)*dx + (0.5)*dx;
              p_solucao[i,j] = (sinh(pi*y)/(sinh(pi)))*sin(pi*x);
              triplets[i,j] = [x, y, p_solucao[i,j]];
#              println("x = ", x)
#              println("y = ", y)
          end
      end

      p_staggered = zeros(n,n); #valor da pressão no meio da célula
      f = zeros(n,n); #valor da f também no meio da célula
      left = zeros(1,n);
      right = zeros(1,n);
      upper = zeros(1,n);
      lower = zeros(1,n);
      for i in -1:n-2
          x = (i+0.5)*dx;
          upper[i+2] = sinpi(x);
      end
      

    solvePoisson!(n, p_staggered, f, left, right, upper, lower)
    pn = zeros(n-2, n-2);
    staggered2not!(p_staggered,pn,n)
#    println(p_solucao - pn)
    return maximum(abs(p_solucao - pn))
end

function staggered2not!(p, pn, n) 
	for i in 2:n-1 
		for j in 2:n-1 
			pn[i-1,j-1] = p[i,j]
		end 
	end
end
