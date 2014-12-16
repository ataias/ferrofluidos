using Gadfly

f = open("errorUVvsDx2.dat");
l = length(13:4:69);

dx2 = zeros(l);
error_u = zeros(l);
error_v = zeros(l);

for i in 1:l
	dx2[i]       = read(f, Float64);
	error_u[i]   = read(f, Float64);
	error_v[i]   = read(f, Float64); 
end

close(f)

plot(x=dx2, y=error_u,
	Scale.y_sqrt, Geom.point, Geom.smooth,
	Guide.xlabel("dx2"), Guide.ylabel("abs(erro)"), 
	Guide.title("Erro da velocidade em x"))

plot(x=dx2, y=error_v,
	Scale.y_sqrt, Geom.point, Geom.smooth,
	Guide.xlabel("dx2"), Guide.ylabel("abs(erro)"), 
	Guide.title("Erro da velocidade em x"))

