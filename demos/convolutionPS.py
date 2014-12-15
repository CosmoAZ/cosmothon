import sys
import math
import matplotlib.pyplot as plt
from scipy.special import jn as J
from numpy import log10
from numpy import linspace
from numpy import logspace
from scipy.integrate import quad
from scipy.interpolate import interp1d

def main(argv):
    print ''
    print '=== Running convolutionPS.py program ====\n'

    print argv

    larray = []
    Pkarray = []

    f = open('l.txt', 'r')
    for line in f:
	larray.append(float(line))
    f.close()

    g = open('Pk.txt', 'r')
    for line in g:
	Pkarray.append(float(line))
    g.close()

    e_plus = []
    e_minus = []

    #theta = linspace(1/60.0*math.pi/180.0, 100.0/60.0*math.pi/180.0, 30)
    theta = logspace(log10(0.01/60.0*math.pi/180.0), log10(100.0/60.0*math.pi/180.0), 300)

    P_l = interp1d(larray, Pkarray)

    f_x = lambda l,t: l/(2*math.pi)*J(0, l*t)*P_l(l)
    f_x_m = lambda l,t: l/(2*math.pi)*J(4, l*t)*P_l(l)

    """for t in theta:

	e_p = 0.0
	e_m = 0.0

	i = 0

	for i in range(len(larray)-1):#larray[0:len(larray)-1]:
	
		e_p += (larray[i]/(2*math.pi) * J(0, larray[i]*t)*Pkarray[i] + larray[i + 1]/(2*math.pi) *  J(0, larray[i + 1]*t) * Pkarray[i+1])/2.0*(larray[i + 1] - larray[i])
		e_m += (larray[i]/(2*math.pi) * J(4, larray[i]*t)*Pkarray[i] + larray[i + 1]/(2*math.pi) *  J(4, larray[i + 1]*t) * Pkarray[i+1])/2.0*(larray[i + 1] - larray[i])


	e_plus.append(e_p)
	e_minus.append(e_m)"""

    for t in theta:

	e_p, errrop = quad(f_x, min(larray),  max(larray), args=(t,), epsabs = 1e-12)
	e_m, errrom = quad(f_x_m, min(larray),  max(larray), args=(t,), epsabs = 1e-12)

	e_plus.append(e_p)
	e_minus.append(e_m)

    i = 0

    eout = open('e_plus.txt', 'w')
    emout = open('e_minus.txt', 'w')

    for i in range(len(theta)):

	theta[i] *= 60.0*180.0/math.pi
	eout.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\n')
        emout.write(str(theta[i]) + '\t' + str(e_minus[i]) + '\n')

    eout.close()
    emout.close()

###########################################

    xstuff = []

    for l in larray:
	for t in theta:
	    xstuff.append(l*math.pi/(60.0*180.0)*t)

    import numpy as np
    import pylab as py
    import scipy.special as sp

    x = np.linspace(0, max(xstuff), 50000000)

    for v in [0, 4]:
    	py.plot(x, sp.jv(v, x))

    py.xlim((0, max(xstuff)))
    py.ylim((-0.5, 1.1))
    #py.legend(('$\mathcal{J}_0(x)$', '$\mathcal{J}_1(x)$', '$\mathcal{J}_2(x)$',
           #'$\mathcal{J}_3(x)$', '$\mathcal{J}_4(x)$', '$\mathcal{J}_5(x)$'),
           #loc = 0)
    py.xlabel('$x$')
    py.ylabel('$\mathcal{J}_n(x)$')
    #py.title('Plots of the first six Bessel Functions')                                
    py.grid(True)
    py.savefig('besseln0to6.png')                                      
    #py.show()
    py.close()

##########################################

    plt.plot(theta, e_minus, 'r-')
    plt.plot(theta, e_plus, 'b-')
    #plt.semilogy()
    plt.semilogx()
    plt.xlabel('theta')
    plt.ylabel('e_plus')
    plt.title('e_plus vs. theta')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.axis([1, max(theta), -5e-5, 3e-4])#min(e_plus), max(e_plus)])
    plt.grid(True)
    plt.savefig('e_plus.png')
    plt.close()

    """plt.plot(theta, e_minus, 'r-')
    #plt.semilogy()
    plt.semilogx()
    #plt.xlabel('theta')
    #plt.ylabel('e_minus')
    #plt.title('e_minus vs. theta')
    #plt.axis([0, max(theta), min(e_minus), max(e_minus)])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True)
    plt.savefig('e_minus.png')
    plt.close()"""


if __name__ == "__main__":
    main(sys.argv[1:])
