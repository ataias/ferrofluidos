using PoissonDirichlet

function testPoissonDirichlet(n = 101)
    #Inicialize as variáveis de interesse
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
