#Navier Stokes validation

from matplotlib import *
from pylab import *
import struct
import sys

#faltando 170 e depois subtrair do resultado de 200
#depois, fazer a mesma coisa daqui com o resultado da parte magnética

y = 100*abs(abs(array([-0.6583694564678838, -0.6598305142767019, -0.6605565482345396, -0.6609702618352404, -0.6612284529070828, -0.6614004357603461])) - 0.6615666703555553)/0.6615666703555553
x = 1/array([70, 90, 110, 130, 150, 170])**2

#use LaTeX, choose nice some looking fonts and tweak some settings
rc('font', family='serif')
rc('font', size=16)
rc('legend', fontsize=16)
rc('legend', numpoints=1)
rc('legend', handlelength=1)
rc('legend', frameon=False)
rc('xtick.major', pad=7)
rc('xtick.minor', pad=7)
rc('text', usetex=True)
rc('text.latex',
             preamble=[r'\usepackage[T1]{fontenc}',
                       r'\usepackage{amsmath}',
                       r'\usepackage{txfonts}',
                       r'\usepackage{textcomp}'])

close('all')

fig = figure()
ax = fig.add_subplot(111)
ax.xaxis.set_major_locator(MaxNLocator(4))

xlabel('$\Delta x^2$')
ylabel('Error (\%)')
title('Vorticity error for middle-point without magnetism')

# axis([-0.00002, 0.00021, -0.000002, 0.00007])
grid(True)

plot(x,y, 'bo', label="Original data")

A = np.vstack([x, ones(len(x))]).T
m, c = linalg.lstsq(A, y)[0]
plot(x, m*x + c, 'r', label='Fitted line')

savefig("validateHydrodynamicsRe40.png", dpi=300)
# show()
