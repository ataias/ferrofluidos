#!/usr/bin/env python3

import glob, os, re, shutil
from vplot import *
from numpy import *
import time as tt
import math

if __name__ == "__main__":

    for filename in glob.glob("*.dat"):
        start = tt.time()
        #Nome do arquivo tem parâmetros chave para código
        parameters = re.findall(r"[-+]?\d*\.\d+|\d+", filename)

        Re = float(parameters[0])
        n = int(parameters[1])
        Pe = float(parameters[2])
        alpha = float(parameters[3])
        Cpm = float(parameters[4])
        t = float(parameters[5])
        fps = int(parameters[6])

        print("Arquivo: ", filename)
        print("Re = ", Re)
        print("n = ", n)
        print("Pe = ", Pe)
        print("alpha = ", alpha)
        print("Cpm = ", Cpm)
        print("t = ", t)
        print("fps = ", fps)

        print()

        #Criar pasta
        directory = os.path.splitext(filename)[0]
        if not os.path.exists(directory):
            os.makedirs(directory)

        #Mover arquivos
        shutil.move(filename, directory + "/" + filename)
        shutil.move(directory + ".txt", directory + "/" + directory + ".txt")

        #Entrar na pasta
        os.chdir(directory)

        #Abrir arquivo
        f = open(filename, 'rb')

        #Parâmetros derivados
        dx = 1.0/n
        def sideTextV(time):
            text(1.21, 1.0, r'Re = ' + str(Re))
            text(1.21, 0.95, r'$\alpha = ' + str(alpha) + '$')
            text(1.21, 0.9, r'Cpm = ' + str(Cpm))
            text(1.21, 0.85, r'$n = ' + str(n) + '$')
            text(1.21, 0.8, r'$\mathrm{Pe} = ' + str(Pe) + '$')
            title(r'$\mathbf{v}$, ' + 't = ' + '{:10.8f}'.format(time))


        def sideTextM(time):
            text(1.21, 1.0, r'Re = ' + str(Re))
            text(1.21, 0.95, r'$\alpha = ' + str(alpha) + '$')
            text(1.21, 0.9, r'Cpm = ' + str(Cpm))
            text(1.21, 0.85, r'$n = ' + str(n) + '$')
            text(1.21, 0.8, r'$\mathrm{Pe} = ' + str(Pe) + '$')
            title(r'$\mathbf{M}$, ' + 't = ' + '{:10.8f}'.format(time))

        def sideTextVorticty():
            xlabel('$t$')
            ylabel('$\omega$')
            title(r'$\omega(0.5,0.5)$ for ' + \
                    'Re = ' + str(Re) + \
                    r', $\alpha = ' + str(alpha) + '$' +\
                    r', Cpm = ' + str(Cpm) +\
                    r', $n = ' + str(n) + '$' +\
                    r', $\mathrm{Pe} = ' + str(Pe) + '$')

        def titleTextMagnetism():
            title(r'$|\mathbf{M}|(0.5,0.5)$ for ' + \
                    'Re = ' + str(Re) + \
                    r', $\alpha = ' + str(alpha) + '$' +\
                    r', Cpm = ' + str(Cpm) +\
                    r', $n = ' + str(n) + '$' +\
                    r', $\mathrm{Pe} = ' + str(Pe) + '$')


        u = zeros((n,n), dtype=float64)
        v = zeros((n,n), dtype=float64)
        p = zeros((n,n), dtype=float64)
        Hx = zeros((n,n), dtype=float64)
        Hy = zeros((n,n), dtype=float64)
        Mx = zeros((n,n), dtype=float64)
        My = zeros((n,n), dtype=float64)
        phi = zeros((n,n), dtype=float64)
        angles = zeros((n,n), dtype=float64)

        #Como há o frame 0, precisa somar 1 ao número de frames
        numberFrames = int(round(fps*t))+1

        #As velocidades estão salvas nas quinas
        x=linspace(0, 1 - dx, n)
        y=linspace(0, 1 - dx, n)
        x, y=meshgrid(x, y)
        time = 0.0

        os.makedirs("png")
        os.chdir("png")

        #Note que o tempo inicial salvo não é 0
        #soma-se 1 em int(t*fps+1) devido ao tempo 0
        tvector = linspace(0, t, int(t*fps+1))
        vortc = zeros(size(tvector))

        modM = zeros(size(tvector))
        phaseM = zeros(size(tvector))
        phaseDiffMH = zeros(size(tvector))

        for i in range(0, numberFrames):
            time = i/fps
            # batchRead(u, v, p, Hx, Hy, Mx, My, phi, angles, n, f)
            #1 = SEEK_CUR -> offset relative to current position

            #Criar gráfico para velocidades
            readMatrix(u, f, n)
            readMatrix(v, f, n)
            plotStreamFrame(u, v, x, y, n, sideTextV, time, "v" + str(i).zfill(4) + ".png")
            print("Criada imagem para campo de velocidades em t = ", time)
            # -------------------------------

            #Evolução da vorticidade no meio
            #Este n é o tamanho da malha não escalonada.
            #No código em Julia, n' = 7 é o tamanho da malha escalonada
            #mas aqui n = n' - 2 = 5, por exemplo
            #Lá, o ponto central é dado por int(n'/2)
            #Aqui, o mesmo é válido, apesar de usar n ao invés de n',
            # pois a indexação começa de 0 em python
            c = int(n/2)
            vortc[i] = ((v[c+1,c]-v[c-1,c]) - (u[c,c+1]-u[c,c-1]) )/(2*dx)
            # -------------------------------


            f.seek(n*n*8*1, 1) #pula pressão p

            readMatrix(Hx, f, n)
            readMatrix(Hy, f, n)

            readMatrix(Mx, f, n)
            readMatrix(My, f, n)
            plotStreamFrame(Mx, My, x, y, n, sideTextM, time, "M" + str(i).zfill(4) + ".png")
            print("Criada imagem para campo M em t = ", time)

            #Evolução do magnetismo no meio
            modM[i] = sqrt((Mx[c,c])**2 + (My[c,c])**2)
            phaseM[i] = math.degrees(math.atan2(My[c,c], Mx[c,c]))
            phaseDiffMH[i] = math.degrees(math.atan2(Hy[c,c], Hx[c,c])) - phaseM[i]

            f.seek(n*n*8*2, 1)

        #Gráfico da vorticidade no meio, evoluindo no tempo
        plotPointEvolution(tvector, vortc, sideTextVorticty, "vort" + directory + ".png")
        print("Criada imagem para evolução temporal da vorticidade em (0.5, 0.5)")

        plotMEvolution(tvector, modM, phaseM, phaseDiffMH, titleTextMagnetism,\
                        "Magnetism" + directory + ".png")
        print("Criada imagem para evolução temporal do magnetismo em (0.5, 0.5)")

        os.chdir("..")
        f.close()

        #Voltar a pasta anterior
        os.chdir("..")
        end = tt.time()
        print("Tempo gasto em " + directory + " foi de ", end - start, " segundos")
