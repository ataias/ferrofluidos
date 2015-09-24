using Transient
using NavierTypes

#Executar
#julia -p 1 simulate.jl &
#-p 1 -> indica que julia deve utilizar somente um processo

t = 6.0;
Re = 50.0
n = 102
divFactor = 1.25
dt = getDt(n,Re, divFactor)

save = true;

#canto esquerdo
a = +0.00;
b = -0.05
step = 1


main = homedir()*"/Documents/simulacao"
#Create working folder
wfolder = main * "/data" #working folder
try rm(wfolder, recursive=true) end
mkdir(wfolder)

function simulation(n, dt, Re, c1, Cpm, alpha, a, b)
  fname = "Re" * string(int(Re)) * "Pe" * string(round(1/c1, 4))
  fname *= "α" * string(round(alpha, 4))
  fname *= "Cpm" * string(int(Cpm))

  txtfilename  = fname * ".txt"
  datafilename = fname * ".dat"


  cd(wfolder)
  #redirect standard output
  f = open(txtfilename, "w")
  redirect_stdout(f)

  #Simulate
  @time transient(n, dt, Re, t, Cpm, alpha, a, b, save, c1, datafilename)
  cd("..")

  close(f)
end #end simulate function

Pe = [0.1, 1, 100]
alpha = [0.1, 100]
Cpm = [1, 10]
for P in Pe
  for α in alpha
    for C in Cpm
      simulation(n, dt, Re, 1/P, C, α, a, b)
    end
  end
end
