module Magnetism

function getForce!(n, Cpm, Hx, Hy, Mx, My, fx, fy)
    Myt = 0.0
    Mxt = 0.0
    
    for i in 2:n-1
        for j in 2:n-1
            Myt = (My[i,j]+My[i,j+1]+My[i-1,j]+My[i-1,j+1])/4
            fx[i,j]  = Cpm*Mx[i,j]*( Hx[i,j+1]-Hx[i,j-1])/(2*dx)
            fx[i,j] += Cpm*Myt[i,j]*(Hx[i+1,j]-Hx[i-1,j])/(2*dx)
            
            Mxt = (Mx[i,j]+Mx[i+1,j]+Mx[i,j-1]+M[i+1,j-1])/4
            fy[i,j]  = Cpm*Mxt[i,j]*(Hy[i+1,j]-Hy[i-1,j])/(2*dx)
            fy[i,j] += Cpm*My[i,j]*( Hy[i,j+1]-Hy[i,j-1])/(2*dx)
        end
    end
    
end

end