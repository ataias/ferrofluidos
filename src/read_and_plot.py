#!/usr/bin/env python3

import h5py
import glob, os, re, shutil
from plot import *
from numpy import *
import time as tt
import math

#Este arquivo cria todas as imagens para os arquivos .h5 na pasta na qual ele for executado

#f is an hdf5 file
def printFileData(info):
    print("Date: ", info['date'].value)
    print("Re = ", info['Re'].value)
    print("n = ", info['n'].value)
    print("Pe = ", info['Pe'].value)
    print("alpha = ", info['alpha'].value)
    print("Cpm = ", info['Cpm'].value)
    print("t = ", info['t'].value)
    print("dt = ", info['dt'].value)
    print("fps = ", info['fps'].value)
    print("Simulation complete: ", bool(info['completed'].value))

    print()

if __name__ == "__main__":

    for filename in glob.glob("*.h5~"):
        start = tt.time()
        #Criar pasta
        directory = os.path.splitext(filename)[0]
        if not os.path.exists(directory):
            os.makedirs(directory)

        #Mover arquivos
        shutil.move(filename, directory + "/" + filename)
        shutil.move(directory + ".txt~", directory + "/" + directory + ".txt")

        #Entrar na pasta
        os.chdir(directory)

        #Abrir arquivo
        # f = open(filename, 'rb')
        f = h5py.File(filename, "r")

        info = f.get("info")
        frames = f.get("frames")
        printFileData(info)

        Re = info['Re'].value
        n = info['n'].value
        Pe = info['Pe'].value
        alpha = info['alpha'].value
        Cpm = info['Cpm'].value
        t = info['t'].value
        fps = info['fps'].value

        #Parâmetros derivados
        dx = 1.0/n

        def setupSideText():
            text(1.21, 1.0, r'Re = ' + str(Re))
            text(1.21, 0.95, r'$\alpha = ' + str(alpha) + '$')
            text(1.21, 0.9, r'Cpm = ' + str(Cpm))
            text(1.21, 0.85, r'$n = ' + str(n) + '$')
            text(1.21, 0.8, r'$\mathrm{Pe} = ' + str(Pe) + '$')

        def sideTextV(time):
            setupSideText()
            title(r'$\mathbf{v}$, ' + 't = ' + '{:10.8f}'.format(time))


        def sideTextM(time):
            setupSideText()
            title(r'$\mathbf{M}$, ' + 't = ' + '{:10.8f}'.format(time))

        def sideTextF(time):
            setupSideText()
            title(r'$\mathbf{F}$, ' + 't = ' + '{:10.8f}'.format(time))

        def sideTextH(time):
            setupSideText()
            title(r'$\mathbf{H}$, ' + 't = ' + '{:10.8f}'.format(time))

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
        fx = zeros((n,n), dtype=float64)
        fy = zeros((n,n), dtype=float64)

        #Como há o frame 0, precisa somar 1 ao número de frames
        numberFrames = frames['lastFrame'].value+1

        #As velocidades estão salvas nas quinas
        x=linspace(0, 1 - dx, n)
        y=linspace(0, 1 - dx, n)
        x, y=meshgrid(x, y)
        time = 0.0

        os.makedirs("png")
        os.chdir("png")

        #Note que o tempo inicial salvo não é 0
        #soma-se 1 em int(t*fps+1) devido ao tempo 0
        tvector = linspace(0, t, numberFrames)
        vortc = zeros(size(tvector))

        modM = zeros(size(tvector))
        phaseM = zeros(size(tvector))
        phaseDiffMH = zeros(size(tvector))

        for i in range(0, numberFrames, 5):
            time = frames['t/' + str(i)].value

            #Criar gráfico para velocidades
            u = frames['velocity/x/' + str(i)].value
            v = frames['velocity/y/' + str(i)].value
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

            Hx = frames['H/x/' + str(i)].value
            Hy = frames['H/y/' + str(i)].value

            pngName = "H" + str(i).zfill(4) + ".png~"
            plotStreamFrame(Hx, Hy, x, y, n, sideTextH, time, pngName)
            print("Criada imagem para campo H em t = ", time)
            shutil.move(pngName, os.path.splitext(pngName)[0] + ".png")

            # readMatrix(Mx, f, n)
            Mx = frames['M/x/' + str(i)].value
            # readMatrix(My, f, n)
            My = frames['M/y/' + str(i)].value

            pngName = "M" + str(i).zfill(4) + ".png~"
            plotStreamFrame(Mx, My, x, y, n, sideTextM, time, pngName)
            print("Criada imagem para campo M em t = ", time)
            shutil.move(pngName, os.path.splitext(pngName)[0] + ".png")

            #Evolução do magnetismo no meio
            modM[i] = sqrt((Mx[c,c])**2 + (My[c,c])**2)
            phaseM[i] = math.degrees(math.atan2(My[c,c], Mx[c,c]))
            phaseDiffMH[i] = math.degrees(math.atan2(Hy[c,c], Hx[c,c])) - phaseM[i]

            # Força
            Fx = frames['F/x/' + str(i)].value
            Fy = frames['F/y/' + str(i)].value
            pngName = "F" + str(i).zfill(4) + ".png~"
            plotVectorFrame(Fx, Fy, x, y, n, sideTextF, time, pngName)
            print("Criada imagem para campo F em t = ", time)
            shutil.move(pngName, os.path.splitext(pngName)[0] + ".png")

        #Gráfico da vorticidade no meio, evoluindo no tempo
        pngName = "vort" + directory + ".png~"
        plotPointEvolution(tvector, vortc, sideTextVorticty, pngName)
        print("Criada imagem para evolução temporal da vorticidade em (0.5, 0.5)")
        shutil.move(pngName, os.path.splitext(pngName)[0] + ".png")

        pngName = "Magnetism" + directory + ".png"
        plotMEvolution(tvector, modM, phaseM, phaseDiffMH, titleTextMagnetism,\
                        pngName)
        print("Criada imagem para evolução temporal do magnetismo em (0.5, 0.5)")
        shutil.move(pngName, os.path.splitext(pngName)[0] + ".png")

        os.chdir("..") #out of png directory
        f.close()

        #Let simulation file be uploaded and then deleted by dropbox daemon
        shutil.move(filename, directory + ".h5")

        #Voltar a pasta anterior
        os.chdir("..") #out of simulation directory and back to "simulations"
        end = tt.time()
        print("Tempo gasto em " + directory + " foi de ", end - start, " segundos")
