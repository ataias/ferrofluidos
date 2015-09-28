#Executar
#julia -p 0 simulate.jl &

# É importante que o número de processos seja definido logo no início
# Caso contrário, nem todos os módulos serão corretamente carregados nos processos
if nprocs() - 1 < CPU_CORES
  addprocs(CPU_CORES - nprocs() + 1)
end

println("nprocs() = ", nprocs())

@everywhere using Transient
@everywhere using NavierTypes

t = 30.0
Re = 50.0
n = 102
divFactor = 1.25
dt = getDt(n,Re, divFactor)

save = true;

#canto esquerdo
a = +0.00;
b = -0.05
step = 1

fps = 2

@everywhere main = homedir()*"/Documents/simulacao"
#Create working folder
@everywhere wfolder = main * "/correct" #working folder
try rm(wfolder, recursive=true) end
mkdir(wfolder)

@everywhere function simulation(n, dt, Re, c1, Cpm, alpha, a, b, t, fps, save)
  # println("pwd() = ", pwd())
  fname = "Re" * string(int(Re)) * "N" * string(int(n-2))
  fname *= "Pe" * string(round(1/c1, 4))
  fname *= "α" * string(round(alpha, 4))
  fname *= "Cpm" * string(int(Cpm))
  fname *= "T" * string(t)
  fname *= "fps" * string(fps)

  txtfilename  = fname * ".txt"
  datafilename = fname * ".dat"


  cd(wfolder)
  #redirect standard output
  f = open(txtfilename, "w")
  redirect_stdout(f)

  #Simulate
  @time transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, datafilename)
  cd("..")

  close(f)
end #end simulate function

Pe = [0.1, 1, 5]
alpha = [0.1, 100]
Cpm = [1, 10]
total = length(Pe) * length(alpha) * length(Cpm)
i = 0
@sync begin
  for P in Pe
    for α in alpha
      for C in Cpm
        @spawnat int(i % CPU_CORES + 2) begin
          simulation(n, dt, Re, 1/P, C, α, a, b, t, fps, save)
        end
        i = i + 1
      end
    end
  end
end
