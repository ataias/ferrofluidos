using NavierStokes
#need to execute this directly... otherwise those variables are all local
n = 6;
p = zeros(n,n);
u = zeros(n,n);
u[:,n] = 2*ones(typeof(1.0),1,n)
println(u)
v = zeros(n,n);
u_old = copy(u);
v_old = copy(v);
fx = zeros(n,n);
fy = zeros(n,n);
mu = 1.0;
rho = 1.0;
Re = 1.0/mu;
# uB = [1.0 1 1 1 1 1];
uB = ones(typeof(1.0),1,n)
dt = getDt(n, Re);


if(!isdXok(1/mu, n))
	println("There is a problem with your dx. Increase n.\n");
	exit(0)
end
@time solve_navier_stokes!(n, dt, mu, rho, p, u, v, u_old, v_old, fx, fy, uB)
println(u)
println(v)