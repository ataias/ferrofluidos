#Plotting Poisson validation

from matplotlib import *
from pylab import *
import struct
import sys

y = 100*array([0.00019758602842993855, 4.958384660265491e-5, 1.2404251665182331e-5, 7.938993296777164e-6, 5.0813968964918965e-6, 3.1017008645523036e-6])/0.006890576458974079 #valor real é o denominador
x = 1/array([40, 80, 160, 200, 250, 320])**2

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

xlabel('$\Delta x^2$')
ylabel('Error (\%)')
title('Error in Poisson solution for middle point')

# axis([-0.00002, 0.00021, -0.000002, 0.00007])
grid(True)

plot(x,y, 'bo', label="Original data")

A = np.vstack([x, ones(len(x))]).T
m, c = linalg.lstsq(A, y)[0]
plot(x, m*x + c, 'r', label='Fitted line')

savefig("validatePoissonP2.png", dpi=300)
# show()
