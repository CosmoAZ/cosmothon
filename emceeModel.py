import sys
import emcee
from numpy import logspace
import shearPSFile
import lensingClassFile
import math
from matplotlib import pyplot as plt
import time
import numpy
from numpy import isfinite
from numpy import *

theta = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 29)
edat = [0.000174, 0.000239, 0.000146, 0.000157, 0.000145, 0.000112, 0.000094, 0.000124, 0.000121, 0.000129, 0.000112, 0.000111, 0.000093, 0.000082, 0.000058, 0.000043, 0.000026, 0.000015, 0.000010, 0.000009, 0.000005, 0.000006, -0.000007, -0.000012, -0.000018, -0.000018, -0.000013, -0.000011, 0.000004]
print edat


lcp = lensingClassFile.lensingClass(0.7, 0.3, 0.7)
F2redshifts = lcp.getRedshiftsFromDLS()
print "Calculating p(z)"
orgList, numList = lcp.organizeRedshifts(F2redshifts)
nlistp = lcp.n_list(numList, orgList, F2redshifts)
print "nlistp: "
print nlistp
print "Finished calculating p(z)"
print "numlist:"
print numList
print "orgList:"
print orgList    	


print nlistp

#sys.exit()

def lnprior(phi):
	omegamat = phi[0]
	omegaDE = phi[1]
	if omegamat < 0 or omegamat > 1 or omegaDE < 0 or omegaDE > 1:
		return -inf
	return 0

def lnlike(phi, theta, edat, F2red, nlistp):

	start = time.time()

	omegamat = phi[0]
	omegaDE = phi[1]

	#print("\nshearClacClass:  Omegamat:\t" + str(omegamat) + "  OmegaDE:\t" + str(omegaDE))
	#e_plus = runshear.runshearPS(p[0], p[1]) 
    	shearPSClass = shearPSFile.shearCalcClass(1, 0, 0, 0.673, omegamat, omegaDE, nlistp)
	Pk_l, l = shearPSClass.shearPSPrePrep(F2red)

    	convolutionClassObj = shearPSFile.convolutionClass(1, 0, 0, 0.673, omegamat, omegaDE)
    	e_plus, e_minus, theta = convolutionClassObj.convolutionCalc(Pk_l, l)

	end = time.time()

	#f = open('stats2.txt', 'w')
	#f.write(str(theta) + '\n')
	#f.close()

	int1 = []
	int2 = []

	e_prime = e_minus
	for t in theta:
		i = 0
		for eta in theta:
			if (eta >= theta):
				int1.append(e_minus[i]/eta)
				int2.append(t*t*e_minus[i]/(eta*eta*eta))
			i += 1


	
	sum=0
	for i in range(len(edat)):
	    sum += pow(edat[i] - e_plus[i], 2.0)
	return math.exp(-sum)


def lnprob(p, theta, edat, F2red, nlistp):

	lp = lnprior(p)
	if not isfinite(lp):
		return -inf
	return lp + lnlike(p, theta, edat, F2red, nlistp)


#import scipy.optimize as opt
#nll = lambda *args: -lnlike(*args)
#result = opt.minimize(nll, [0.3, 0.7],
#                     args=(theta, edat, F2redshifts))
#print result['x']

trueomegamat= 0.315
trueomegaDE = 0.685
trueh = 0.673

shearPSClass = shearPSFile.shearCalcClass(1, 0, 0, trueh, trueomegamat, trueomegaDE, nlistp)
Pk_l, l = shearPSClass.shearPSPrePrep(F2redshifts)

convolutionClassObj = shearPSFile.convolutionClass(1, 0, 0, trueh, trueomegamat, trueomegaDE)
e_plus, e_minus, theta = convolutionClassObj.convolutionCalc(Pk_l, l)
print e_plus
sys.exit()

Nwalkers, Ndim = 4, 2
pos = [[trueomegamat, trueomegaDE] for i in range(Nwalkers)]#+ 1e-1*numpy.random.randn(Ndim) 

#random.seed(time.time())

#p0 = [0.3, 0.7]# + random.randn(Ndim) for i in range(Nwalker)]

sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnprob, args=(theta, edat, F2redshifts, nlistp), threads=2)

pos,prob,state = sampler.run_mcmc(pos, 1)

samples = sampler.chain[:, :, :].reshape((-1, Ndim))

print pos
print prob
print state
print sampler


res=plt.plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#plt.axhline(alpha_true, color='blue')
plt.axis([0, 25, 0, 1])
plt.savefig("plt.png")
plt.close()

import triangle
#tmp = triangle.corner(sampler.flatchain, labels=['omegamat','omegaDE'], truths=[0.25, 0.75])

fig = triangle.corner(samples)#, labels=["$m$", "$b$"],truths=[0.277, 0.723])
fig.savefig("triangle.png")

#res=plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#axhline(alpha_true, color='blue')


