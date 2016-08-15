using ArgParse
using Transient
using NavierTypes

#Exemplo:
# julia frames.jl --help para mostrar ajudar
# 		julia frames.jl 50 10.0 2.5 0.5 0.8 0.8 0.0 -0.05

#Antes, era necessário entrar com "52" para um tamanho de malha de 50, por conta de 2 colunas/linhas extras da malha escalonada
#Modifiquei isso, mas sem modificar todo o resto do código então precisa de compensação neste
EXTRA_STAGGERED_GRID = 2

function parse_commandline()
    s = ArgParseSettings(description = "Este programa simula as equações de Navier Stokes em uma cavidade sob ação de um campo magnético aplicado.",
                         version = "1.0",
                        add_version = true)

    @add_arg_table s begin
        "--divFactor", "-d"
          help = "Fator de divisão para obter dt, deve ser maior que 1.0"
          arg_type = Float64
          default = 1.25

        "--fps", "-f"
          help = "Número de frames por segundo que devem ser salvos"
          arg_type = Int
          default = 30

        "--save"
          help = "Salvar dados da simulação em arquivo HDF5"
          action = :store_true

        "n"
          help = "Tamanho da malha, deve ser maior que o número de Reynolds da simulação"
          arg_type = Int
          required = true

        "Re"
          help = "Número de Reynolds"
          arg_type = Float64
          required = true

        "t"
          help = "Tempo total de simulação, deve ser maior que 0. A unidade é adimensional."
          arg_type = Float64
          required = true

        "c1"
          help = "Constante relacionada a 1/tau adimensionalizado"
          arg_type = Float64
          required = true

        "Cpm"
          help = "Variável adimensional que influencia diretamente na magnitude da força magnética."
          arg_type = Float64
          required = true

        "alpha"
          help = "Constante relacionada ao campo magnético inicial obtido com a função de Langevin. É relacionado com a intensidade do campo magnético."
          arg_type = Float64
          required = true

        "a"
          help = "Posição x do centro do campo magnético aplicado, deve ser menor que 0 ou maior que 1 para estar fora do quadrado unitário"
          arg_type = Float64
          default = -0.2

        "b"
          help = "Posição y do centro do campo magnético aplicado, deve ser menor que 0 ou maior que 1 para estar fora do quadrado unitário"
          arg_type = Float64
          default = -0.2

    end

    return parse_args(s)
end

function assertParsed(parsed_args)
  @assert(parsed_args["n"] > 10 && parsed_args["n"] > parsed_args["Re"], "n: Use um tamanho de malha maior.")
  @assert(parsed_args["t"] > 0, "t: Tempo de simulação deve ser maior que zero.")
  @assert(parsed_args["Re"] > 0, "Re deve ser maior que zero.")
  @assert(parsed_args["divFactor"] > 1, "Fator de divisão deve ser maior que 1.")
  @assert(!((0 < parsed_args["b"] < 1) && (0 < parsed_args["a"] < 1)), "Centro do campo magnético deve ser definido fora da cavidade.")
  @assert(parsed_args["fps"] > 0, "fps deve ser maior que 0. Caso não queira salvar, não configure ative a opção --save")
  #TODO criar asserts para c1, Cpm, alpha
end

function main()

  #Command line parsing
  parsed_args = parse_commandline()
  # assertParsed(parsed_args)

  n = parsed_args["n"] + EXTRA_STAGGERED_GRID

  #Name of files that will be saved
  basefilename = "Re" * string(round(Int,parsed_args["Re"])) * "N" * string(n-2)
  filename = basefilename * ".txt~"
  datafilename = basefilename *".h5~"

  file = open(filename, "w")
  redirect_stdout(file)

  #TODO: Add `file` as a parameter to transient and flush it more frequently so that the user can see the evolution of the problem in the .txt~ file.

  @time transient(n,
                  getDt(n, parsed_args["Re"], parsed_args["divFactor"]),
                  parsed_args["Re"],
                  parsed_args["t"],
                  parsed_args["Cpm"],
                  parsed_args["alpha"],
                  parsed_args["a"],
                  parsed_args["b"],
                  parsed_args["save"],
                  parsed_args["c1"],
                  parsed_args["fps"],
                  datafilename,
                  1.0)
  close(file)
end

main()
