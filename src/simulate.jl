using Transient
using NavierTypes

#Executar
#julia -p 1 simulate.jl 13 &
#-p 1 -> indica que julia deve utilizar somente um processo

t = 3.0;
Re = vcat(1,10:10:100)
#n = map(x -> (x + 10) - (x+10) % 10, Re + 20) + 2 #se quiser ter malhas diferentes...
n = ones(Int64, size(Re))*110 + 2 #Considerando que a simulação ocorrerá até Re 100, uma malha de 110x110 deve ser ok para todos
dt = zeros(size(n))
divFactor = 1.25
for i=1:size(n,1)
    dt[i] = getDt(n[i],Re[i], divFactor)
end
chi = 0.5;
Cpm = 0.8;
save = true;
#canto esquerdo
a = +0.00;
b = -0.05
step = 1


main = homedir()*"/Documents/simulacao"

function simulation(n, dt, Re, chi, Cpm, gamma, a, b, mag)
  folder = string(mag) * "Re" * string(int(Re)) * "N" * string(n-2)
  filename = folder * ".txt"
    
  #Create working folder
  wfolder = main * "/" * folder #working folder
  try rm(wfolder, recursive=true) end
  mkdir(wfolder)
  cd(wfolder)
  
  #redirect standard output
  f = open(filename, "w")
  redirect_stdout(f)
    
  #Simulate
  @time transient(n, dt, Re, t, chi, Cpm, gamma, a, b, save)
  
  #File operations
  mkdir("png")
  cd("./png")
  simulfile =  "N"*string(n - 2)*".dat"
  mv("../"*simulfile, "./"*simulfile)
  cd("..")
  println("\n Execute o script python de dentro da pasta png")
  close(f)
end

s = int(ARGS[1]) #simulation number

#Há 11 simulações de cada tipo
#Argumentos entre 1 e 11 realizam as simulações nomag
#Argumentos entre 12 e 22 realizam as simulações mag
if 1 <= s <= size(Re,1)
    simulation(n[s], dt[s], Re[s]*1.0, 0.0, 0.0, 0.0, a, b, "n") #no magnetism
elseif s <= int(2*size(Re,1))
    simulation(n[s%11], dt[s%11], Re[s%11]*1.0, 0.5, 0.8, 3.5, a, b, "w") #with magnetism
else
    println("Invalid simulation number")
end