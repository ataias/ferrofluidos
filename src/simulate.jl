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
a = -1.0
b = -1.0

fps = 2

@everywhere main = homedir()*"/Documents/simulacao"
#Create working folder
@everywhere wfolder = main * "/superparamagnetic" #working folder
try rm(wfolder, recursive=true) end
mkdir(wfolder)

@everywhere function simulation(n, dt, Re, Cpm, chi, a, b, t, fps, save)
  # println("pwd() = ", pwd())
  fname = "Re" * string(int(Re)) * "N" * string(int(n-2))
  fname *= "χ" * string(round(alpha, 4))
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
  @time transient(n, dt, Re, t, Cpm, chi, a, b, save, fps, datafilename)
  cd("..")

  close(f)
end #end simulate function

chi = [0.1, 100]
Cpm = [1, 10]
i = 0

@sync begin
  for χ in chi
    for C in Cpm
      @spawnat int(i % CPU_CORES + 2) begin
        simulation(n, dt, Re, C, χ, a, b, t, fps, save)
      end
      i = i + 1
    end
  end
end
