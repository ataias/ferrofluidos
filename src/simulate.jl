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
n = 102
divFactor = 1.25

save = true;

#canto esquerdo
a = -0.2
b = -0.2

fps = 2

@everywhere main = homedir()*"/Documents/simulacao"
#Create working folder
@everywhere wfolder = main * "/shliomis" #working folder
try rm(wfolder, recursive=true) end
mkdir(wfolder)

@everywhere function simulation(n, dt, Re, c1, Cpm, alpha, a, b, t, fps, save, convective=false)
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
  convec = 0
  if convective
    convec = 1.0
  end
  @time transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, fps, datafilename, convec)
  cd("..")

  close(f)
end #end simulate function

Re = [1.0, 10.0, 50.0]
Pe = [0.1, 5.0]
alpha = [0.1, 50.0]
Cpm = [1.0, 10.0]

i = 0

@sync begin
  #Simulações para casos não magnéticos
  for R in Re
    dt = getDt(n, R, divFactor)
    @spawnat int(i % CPU_CORES + 2) begin
      simulation(n, dt, R, 1, 0, 0, a, b, 10.0, fps, save, false)
    end #end spawnat
    i = i + 1
 end #end for R in Re

 #Simulações para casos magnéticos sem termo convectivo
  for R in Re
    for P in Pe
      for α in alpha
        for C in Cpm
          dt = getDt(n, R, divFactor)
          @spawnat int(i % CPU_CORES + 2) begin
            simulation(n, dt, R, 1/P, C, α, a, b, t, fps, save, false)
            simulation(n, dt, R, 1/P, C, α, a, b, t, fps, save, true)
          end #end spawnat
          i = i + 1
        end #end for C in Cpm
      end #end for α in alpha
    end #end for P in Pe
  end #end for R in Re
end #end @sync begin
