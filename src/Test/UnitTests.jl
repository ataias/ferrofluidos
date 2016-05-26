using Transient
using NavierTypes
using Gadfly
#N is a vector of all Ns that are going to be used!
# Run everything up to 1 or 2 times
function NavierStokesVaryingNTest(N=[52])
  Re = 10
  c1 = 1 # para caso não magnético, basta não ser 0, já que dá problema na divisão
  Cpm = 0 # forçar caso não magnético
  alpha = 1.0 # semelhantemente
  a, b = 0.1, 0.1 # precisa?
  t = 1.0 # é bom variar também
  fps = 50
  divFactor = 1.25
  save = false #testes unitários só se interessam no ponto do meio
  datafilename = "null"

  results = Dict{Int,Any}()

  for n in N
    dt = getDt(n, Re, divFactor)
    @time results[n] = transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, datafilename)
  end

  # Need to know length of an array in results (depends on fps)
  points = length(results[N[1]][:t])
  midTPoint = round(Int, points/2)

  t = zeros(length(N))
  u = zeros(length(N))
  v = zeros(length(N))
  vortc = zeros(length(N))

  for i in range(1,length(N))
    k = N[i]
    t[i] = results[k][:t][midTPoint]
    println("Time of N = " * string(N[i]) * " is " * string(t[i]))
    u[i] = results[k][:u][midTPoint]
    v[i] = results[k][:v][midTPoint]
    vortc[i] = results[k][:vortc][midTPoint]
  end

  dx = [(1/(i-2))^2 for i in N]
  base_u = u[length(N)]
  du = [abs(x - base_u) for x in u]
  plotu = plot(layer(x=dx, y=du, Geom.point), Guide.XLabel("dx2"), Guide.YLabel("Error"), Guide.Title("Evolution of error according to mesh size step"))
  draw(PDF("plotu.pdf", 8inch, 6inch), plotu)
end

function main()
  NavierStokesVaryingNTest(2 + [50, 60, 70, 80, 90, 100])
end

main()
