using Transient
using NavierTypes

t = 3.0;
Re = vcat(1,10:10:30)
n = map(x -> (x + 10) - (x+10) % 10, Re + 20) + 2
dt = zeros(size(n))
divFactor = 1.25
for i=1:size(n,1), j=n, k=Re
    dt[i] = getDt(j,k, divFactor)
end
chi = 0.5;
Cpm = 0.8;
save = true;
#canto esquerdo
a = +0.00;
b = -0.05
step = 1


main = homedir()*"/Documents/simulacao"

function simulation(j, k, i, chi, Cpm, gamma, a, b, mag)
  folder = string(mag) * "Re" * string(i)
  filename = folder * ".txt"
  f = open(filename, "w")
  redirect_stdout(f)
  @time transient(j, k, i, t, chi, Cpm, gamma, a, b, save)
    #File operations
  try rm(main * folder, recursive=true) end
  try mkdir(main * "/" * folder) end
  try rm("png", recursive=true) end
  try mkdir("png") end
  cd("./png")
  simulfile =  "N"*string(j - 2)*".dat"
  mv("../" * simulfile, "./" * simulfile)
  run(`../vectorField.py $((j-2)) $t $step $chi $Cpm $i $gamma`)
  mv("simulfile", main*folder)
  cd("..")
  close(f)
  mv(filename, main*folder)
  mv("png", main*folder*"png")
end

for i=Re, j=n, k=dt
  simulation(j, k, i*1.0, 0.0, 0.0, 0.0, a, b, "nomag") #without magnetism
end

for i=Re, j=n, k=dt
  simulation(j, k, i, 0.5, 0.8, 3.5, a, b, "mag") #with magnetism
end