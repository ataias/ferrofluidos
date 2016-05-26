using Transient
using NavierTypes
using Gadfly
@everywhere using Transient
@everywhere using NavierTypes
@everywhere using Gadfly

#N is a vector of all Ns that are going to be used!
# Run everything up to 1 or 2 times
function NavierStokesVaryingNTest(;Re=10, divFactor=1.25, t=1.0, N=[52])

  # Forcing non-magnetic case
  c1 = 1
  Cpm = 0
  alpha = 1.0
  a, b = 0.1, 0.1 # precisa?
  # --------------------------

  # Frames per second
  fps = 50

  # In the unit tests, we are not saving the results
  save = false
  datafilename = "null"

  # Dictionary for which key is the mesh size n and its values are dictionaries
  # returne by transient()
  results = Dict{Int,Any}()
  
  println("Iniciando simulação...")
  #Simulações para casos não magnéticos

  @sync begin
    for n in N
      println(n)
      @async begin
        dt = getDt(n, Re, divFactor)
        call = @spawn transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, datafilename, should_print=false)
        println("Tasks n = $(n) started executing")
        results[n] = fetch(call)
        println("Task n = $(n) ended")
      end
    end
  end

  println("Ended executing tasks")
  # Need to know length of an array in results (depends on fps)
  points = length(results[N[1]][:t])
  midTPoint = round(Int, points/2)

  t = zeros(length(N))
  u = zeros(length(N))
  v = zeros(length(N))
  vortc = zeros(length(N))

  for i in range(1,length(N))
    k = N[i]
    # If you want to check whether the time point is actually the same,
    # uncomment the following two lines
    # t[i] = results[k][:t][midTPoint]
    # println("Time of N = " * string(N[i]) * " is " * string(t[i]))

    u[i] = results[k][:u][midTPoint]
    v[i] = results[k][:v][midTPoint]
    vortc[i] = results[k][:vortc][midTPoint]
  end

  dx = [1/(i-2) for i in N]
  dx2 = [x*x for x in dx]
  base_u = u[length(N)]
  du = [abs(x - base_u) for x in u]
  test_plot = plot(layer(x=dx2, y=du, Geom.point), Guide.XLabel("dx2"), Guide.YLabel("Error"), Guide.Title("Evolution of error: Re=" * string(Re) * ", divFactor=" * string(divFactor)))
  draw(PDF("NavierStokesVaryingNTest_Re" * string(round(Int, Re)) * "_divFactor" * string(round(Int,divFactor*100)) * ".pdf", 8inch, 6inch), test_plot)
  println("Finished plotting for test No Magnetism with Re=$(round(Int, Re)) and divFactor=$(round(Int,divFactor*100))")

  # Linear regression of type du = w1 + w2*dx + w3*dx^2
  # We want a and b to be basically 0, this is how we assess this unit test passed
  # How great an epsilon should we choose? I don't know.

  H = [ones(length(dx)) dx dx2]
  y = du
  w = (H'*H)\H'*y

  println("Coefficients are:")
  println("w[1] = $(w[1])")
  println("w[2] = $(w[2])")
  println("w[3] = $(w[3])")
end

function main()
  NavierStokesVaryingNTest(Re=10.0, divFactor=1.25, t=1.0, N=2 + [50, 60, 70])
end

main()
