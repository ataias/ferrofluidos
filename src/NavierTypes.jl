module NavierTypes

export NSEquation, getDt, createNSObject
export VF

using Poisson

#Vector Field
type VF
    x::Array{Float64, 2}
    y::Array{Float64, 2}
end

type NSParams
    n::Int64
    dt::Float64
    Re::Float64
    dx::Float64
    dx2::Float64
    dtx::Float64
    rhodtdx::Float64
end

type NSSystem
    A::SparseMatrixCSC{Float64, Int64}
    chol #CholmodFactor{Float64,Int64}
    b::Array{Float64, 1}
end

#This type holds essential information for the NavierStokes methods
type NSEquation
    params::NSParams
    msystem::NSSystem
    p::Array{Float64, 2}
    v::VF
    v_old::VF
    f::VF
    uB::Array{Float64, 1}
end

function createNSObject(n::Int64, Re::Float64, divFactor::Float64 = 1.25)
    #--------- NSParams ------------
    dx = 1/(n-2)
    if(!isdXok(Re, n))
		println("There is a problem with your dx. Increase n.\n")
		return(1)
	end
    dt = getDt(n, Re, divFactor)
    dx2 = dx*dx
    dtx = dt/(2*dx)
    dtdx = 1/dt/dx
    params = NSParams(n, dt, Re, dx, dx2, dtx, dtdx)

    #-------NSSystem ---------------
    spn = (n-2)*(n-2)
    A = getANeumannSparse(n)
    chol = cholfact(A)
    b = zeros(spn)
    msystem = NSSystem(A, chol, b)

    #------- NSEquation
    p = zeros(n,n)
    v = VF(zeros(n,n), zeros(n,n))
    v_old = VF(zeros(n,n), zeros(n,n))
    f = VF(zeros(n,n), zeros(n,n))
    uB = zeros(n)
    NSObject = NSEquation(params, msystem, p, v, v_old, f, uB)
    NSObject.v.x[1,:] = NaN
    NSObject.v.y[:,1] = NaN
    NSObject.v_old.x[1,:] = NaN
    NSObject.v_old.y[:,1] = NaN

    return NSObject
end

#isdXok
#Esta função analisa se o número de pontos escolhido
#satisfaz a condição de camada limite hidrodinâmica
function isdXok(Re, n)
	dx = 1/(n-2);
	if dx < (1/Re)
		return true
	else
		return false
	end
end

function getDt(n, Re, divFactor=1.25)
	dx = 1/(n-2) #n é o tamanho da malha escalonada
	return min(0.25*Re*dx*dx,dx)/divFactor
end

end #NavierTypes
