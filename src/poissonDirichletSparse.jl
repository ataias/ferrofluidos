function poissonDirichletSparseSolver(n, p, f, A, K, left, right, upper, lower)
    spn = (n-2)*(n-2) #sparse n
    dx = 1/(n-1)
    dx2 = dx*dx

    #we need to solve Ax=b
#    A = getA(n)
#    K = getK(n)
    
    aux = 0.0
#   pCon = zeros(n,n) #conditions are arranged in this matrix
    p[1,:] = left
    p[n,:] = right
    p[:,1] = lower
    p[:,n] = upper
    b = zeros(spn) #remember, we need to solve Ax = b
    
    for i in 1:spn
        aux = 0.0
	    if K[i][1] == 2
	    	aux+=p[1,K[i][2]]
        end
        
	    if K[i][1]+1 == n
	    	aux+=p[n,K[i][2]]
        end
        
	    if K[i][2] == 2
	    	aux+=p[K[i][1],1]
        end
        
	    if K[i][2]+1 == n
	   	    aux+=p[K[i][1],n]
        end
        
	    b[i] = dx2*f[K[i][1], K[i][2]]-aux;
	end
    
    #A was modified in order to be positive definite,
    #before it was negative definite
    x = A\b;
    #arrange result in matrix
    for i in 1:spn
		p[K[i][1],K[i][2]] = -x[i]
	end
end

function getA(n)
    spn = (n-2)*(n-2) #sparse n
    #Definir matriz A
    A = 4*speye(spn, spn)
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
    
    return A
end

function getK(n)
    spn = (n-2)*(n-2)
    K = fill((0,0), (spn,1)) #vetor com os índices dos pontos internos da matriz p
    for i in 2:n-1
	    for j in 2:n-1
	    	K[j+(i-2)*(n-2)-1] = (j,i)
        end
	end
    
    return K
end

function test(n = 101)
    p = zeros(n,n)
    f = zeros(n,n); #valor da f também no meio da célula
    left = zeros(1,n);
    right = zeros(1,n);
    upper = zeros(1,n);
    lower = zeros(1,n);
    dx = 1/(n-1)
    for i in 1:n
        x = (i-1)*dx;
        upper[i] = sinpi(x);
    end
    
    A = getA(n)
    K = getK(n)
    poissonDirichletSparseSolver(n, p, f, A, K, left, right, upper, lower)
    
    p_solucao = zeros(n,n);
    for i in 1:n
        for j in 1:n
            #pressão avaliada no centro da célula
            x = (i-1)*dx;
            y = (j-1)*dx;
            p_solucao[i,j] = (sinh(pi*y)/(sinh(pi)))*sinpi(x);
        end
    end
    
    println("Error = ", maximum(abs(p_solucao-p)))
end
