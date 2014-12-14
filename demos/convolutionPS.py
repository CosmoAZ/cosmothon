import sys
import math
import matplotlib.pyplot as plt
from scipy.special import jv as J
from numpy import log10
from numpy import linspace
from numpy import logspace

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
    theta = logspace(log10(1/60.0*math.pi/180.0), log10(100.0/60.0*math.pi/180.0), 300)

    for t in theta:

	e_p = 0.0
	e_m = 0.0

	i = 0

	for l in larray[0:len(larray)-1]:
	
		e_p += (larray[i]/(2*math.pi) * J(0, larray[i]*t)*Pkarray[i] + larray[i + 1]/(2*math.pi) *  J(0, larray[i + 1]*t) * Pkarray[i+1])/2.0*(larray[i + 1] - larray[i])
		e_m += (larray[i]/(2*math.pi) * J(4, larray[i]*t)*Pkarray[i] + larray[i + 1]/(2*math.pi) *  J(4, larray[i + 1]*t) * Pkarray[i+1])/2.0*(larray[i + 1] - larray[i])

		i += 1

	e_plus.append(e_p)
	e_minus.append(e_m)

    i = 0

    eout = open('e_plus.txt', 'w')
    emout = open('e_minus.txt', 'w')

    for d in theta:

	theta[i] *= 60.0*180.0/math.pi
	eout.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\n')
        emout.write(str(theta[i]) + '\t' + str(e_minus[i]) + '\n')
        i += 1

    eout.close()
    emout.close()

    plt.plot(theta, e_plus, 'b-')
    #plt.semilogy()
    plt.semilogx()
    plt.xlabel('theta')
    plt.ylabel('e_plus')
    plt.title('e_plus vs. theta')
    plt.axis([0, max(theta), min(e_plus), max(e_plus)])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True)
    plt.savefig('e_plus.png')
    #plt.close()

    plt.plot(theta, e_minus, 'r-')
    #plt.semilogy()
    plt.semilogx()
    #plt.xlabel('theta')
    #plt.ylabel('e_minus')
    #plt.title('e_minus vs. theta')
    #plt.axis([0, max(theta), min(e_minus), max(e_minus)])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True)
    plt.savefig('e_minus.png')
    plt.close()


if __name__ == "__main__":
    main(sys.argv[1:])
