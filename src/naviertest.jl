using NavierStokes
#need to execute this directly... otherwise those variables are all local
n = 43;
dx = 1/(n-2)
p = zeros(n,n);
u = zeros(n,n);
uB = zeros(1,n);
for i in -1:n-2
	uB[i+2] = (sinpi(i*dx))^2
	#uB[i+2] = 1
end
u[:,n] = 2*uB
u[1,:] = NaN
	
v = zeros(n,n);
v[:,1] = NaN

u_old = copy(u);
v_old = copy(v);
fx = zeros(n,n);
fy = zeros(n,n);
mu = 1;
rho = 1.0;
Re = 1.0/mu;
dt = getDt(n, Re, 1.5);


if(!isdXok(1/mu, n))
	println("There is a problem with your dx. Increase n.\n");
	exit(0)
end

solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)

# @time for i in 1:2500
# 	solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)
# 	u, u_old = u_old, u
# 	u[:,:] = u_old[:,:]
# end

un = zeros(n-2,n-2)
vn = zeros(n-2,n-2)
pn = zeros(n-2,n-2)
staggered2not!(u, v, p, un, vn, pn, n)