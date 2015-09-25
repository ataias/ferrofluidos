#!/usr/bin/env python3

import glob, os, re, shutil
from streamlines import *

if __name__ == "__main__":

    for filename in glob.glob("*.dat"):
        #Nome do arquivo tem parâmetros chave para código
        parameters = re.findall(r"[-+]?\d*\.\d+|\d+", filename)

        Re = float(parameters[0])
        n = int(parameters[1])
        Pe = float(parameters[2])
        alpha = float(parameters[3])
        Cpm = float(parameters[4])
        t = float(parameters[5])

        print("Arquivo: ", filename)
        print("Re = ", Re)
        print("n = ", n)
        print("Pe = ", Pe)
        print("alpha = ", alpha)
        print("Cpm = ", Cpm)
        print("t = ", t)
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

        u = zeros((n,n), dtype=float64)
        v = zeros((n,n), dtype=float64)
        p = zeros((n,n), dtype=float64)
        Hx = zeros((n,n), dtype=float64)
        Hy = zeros((n,n), dtype=float64)
        Mx = zeros((n,n), dtype=float64)
        My = zeros((n,n), dtype=float64)
        phi = zeros((n,n), dtype=float64)
        angles = zeros((n,n), dtype=float64)

        numberFrames = int(round(180*t))

        x=linspace(0, 1 - dx, n) + (dx/2)
        y=linspace(0, 1 - dx, n) + (dx/2)
        x, y=meshgrid(x, y)
        time = 0.0

        os.makedirs("png")
        os.chdir("png")

        for i in range(0, numberFrames):
            time = (i)/180.0
            try:
                # batchRead(u, v, p, Hx, Hy, Mx, My, phi, angles, n, f)
                #1 = SEEK_CUR -> offset relative to current position

                readMatrix(u, f, n)
                readMatrix(v, f, n)
                plotStreamFrame(u, v, x, y, n, sideTextV, time, "v" + str(i).zfill(4) + ".png")

                f.seek(n*n*8*3, 1)

                readMatrix(Mx, f, n)
                readMatrix(My, f, n)
                plotStreamFrame(Mx, My, x, y, n, sideTextM, time, "M" + str(i).zfill(4) + ".png")

                f.seek(n*n*8*2, 1)


            except:
                print("Exception occurred!")
                break
        os.chdir("..")
        f.close()

        #Voltar a pasta anterior
        os.chdir("..")
