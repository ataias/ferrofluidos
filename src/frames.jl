using NavierStokes
using NavierTypes
using Magnetism
using Poisson

#modo de usar
#n   - ARGS[1] é o tamanho da matriz escalonada
#t   - ARGS[2] é o tempo de simulação, em segundos
#Re  - ARGS[3] é o número de Reynolds
#divFactor  - ARGS[4] dividir o passo de tempo, deve ser maior do que 1
#chi - ARGS[5]
#Cpm - ARGS[6]
#save - ARGS[7] se for 1, salva um arquivo com as matrizes evoluindo no tempo
#a - ARGS[8] é o deslocamento em x da posição central do campo magnético
#b - ARGS[9] é o deslocamento em y da posição central do campo magnético
#gamma - ARGS[10] é a intensidade do campo magnético
#Exemplo:
# 		julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 0.0 3
# Se quiser salvar num arquivo a saída do terminal:
#       julia frames.jl 52 2.5 10.0 1.25 0.5 0.8 0 0.0 0.0 3 > out.txt

n = int(ARGS[1]);
t = float(ARGS[2]);
Re = float(ARGS[3]);
divFactor = float(ARGS[4]);
dt = getDt(n, Re, float(ARGS[4]));
chi = float(ARGS[5]);
Cpm = float(ARGS[6]);
save = bool(int(ARGS[7]));
a = float(ARGS[8]);
b = float(ARGS[9]);
gamma = float(ARGS[10]);

println("Dados sobre simulação:\n n\t= ", n, "\n dx\t= ", 1/(n-2), "\n t\t= ", t, "\n Re\t= ", Re, "\n dt\t= ", dt, "\n ", strftime(time()), "\n");

function factor(i) 
    if i <= 100
        return i/100.0
    else
        return 1.0
    end
end 

#steadyState
#retorna o valor em regime permanente do ponto 0.5, 0.5
#resolve equações para um dado n e Re
#t is time, in seconds, of physical simulation
function steadyState(n, dt, Re, t, chi, Cpm, gamma, save)

	dx = 1/(n-2)
	steps = integer(t/dt)
	c = integer(n/2); #center, non-staggered grid

    NS = createNSObject(n, Re, divFactor)

	for i in -1:n-2
		NS.uB[i+2] = 0*(sinpi(i*dx))^2
	end

    #Condições de contorno, basta editar em v_old
	NS.v_old.x[:,n] = 2*NS.uB

    #Força
	#fx = zeros(n,n)
	#fy = zeros(n,n)

    if(save)
	   file = open("N" * string(n-2) * ".dat", "w")
    end

	un = zeros(n-2,n-2)+1e-15 #malha não escalonada
	vn = zeros(n-2,n-2)+1e-15
	pn = zeros(n-2,n-2)+1e-15

	numberFrames = integer(180*t)
	timeToSave = integer(steps/numberFrames)
    
    #Variáveis para a parte magnética
    phi = zeros(n,n);
    Mx = zeros(n,n);
    My = zeros(n,n);
    left = zeros(n);
    right = zeros(n);
    upper = zeros(n);
    lower = zeros(n);
    for i in -1:n-2
      x = (i+0.5)*dx;
#      upper[i+2] = 1
#      left[i+2] = x
#      right[i+2] = x
       upper[i+2] = sinpi(x)^2
    end
    Hx = zeros(n,n)
    Hy = zeros(n,n)
    
    #non-staggered forms
    Hxn = zeros(n-2, n-2)+1e-15
    Hyn = zeros(n-2, n-2)+1e-15
    phin = zeros(n-2, n-2)+1e-15
    A = getANeumannSparse(n);
#    gamma = 5;
    fHx = (x,y) ->  gamma/(2*pi)*(y-b)/((x-a)^2+(y-b)^2)
    fHy = (x,y) -> -gamma/(2*pi)*(x-a)/((x-a)^2+(y-b)^2)
    
    fact = 0
    if(save)
        write(file, un); write(file, vn);
        write(file, pn);
        write(file, Hxn); write(file, Hyn);
        write(file, phin);
    end
    
	for i in 1:steps
        fact = factor(i)
        for i in -1:n-2
		  NS.uB[i+2] = fact*(sinpi(i*dx))^2
	    end
        getPhi!(n, phi, Mx, My, fHx, fHy, A)
        getMH!(n, chi*fact, phi, Mx, My, Hx, Hy)
        getForce!(n, Cpm, Hx, Hy, Mx, My, NS.f.x, NS.f.y);
        
		solve_navier_stokes!(NS)

		if (i % timeToSave == 0) || (i == steps)
			staggered2not!(NS.v.x, NS.v.y, NS.p, un,  vn,  pn,   n)
            staggered2not!(Hx,     Hy,     phi,  Hxn, Hyn, phin, n)
            if(save)
                write(file, un); write(file, vn);
                write(file, pn);
                write(file, Hxn); write(file, Hyn);
                write(file, phin);
            end
			println("t = ", i*dt)
			println("  u[0.5,0.5]\t= ", un[c,c])
			println("  v[0.5,0.5]\t= ", vn[c,c])
			vortc =  ((vn[c+1,c]-vn[c-1,c]) - (un[c,c+1]-un[c,c-1]))/(2*dx)
			println("  ω [0.5,0.5]\t= ", vortc)
            println("  Pressure values in range\t [", minimum(pn), ", ", maximum(pn), "]")
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

@time steadyState(n, dt, Re, t, chi, Cpm, gamma, save)
