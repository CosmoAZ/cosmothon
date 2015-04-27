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
from pymc.Matplot import plot as samplePlot
import pymc
from pymc import *

theta = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 29)
edat = [0.000174, 0.000239, 0.000146, 0.000157, 0.000145, 0.000112, 0.000094, 0.000124, 0.000121, 0.000129, 0.000112, 0.000111, 0.000093, 0.000082, 0.000058, 0.000043, 0.000026, 0.000015, 0.000010, 0.000009, 0.000005, 0.000006, -0.000007, -0.000012, -0.000018, -0.000018, -0.000013, -0.000011, 0.000004]
print edat
e_low = [5.54E-05, 2.88E-05, 1.63E-05, 2.81E-06, 1.33E-05, 1.84E-05, 2.13E-05, 4.72E-06, 9.09E-06, 1.55E-05, 1.27E-05, 5.85E-06, 8.85E-06, 1.26E-05, 5.29E-06, 7.18E-06, 4.48E-06, 7.37E-06, 5.60E-06, 5.71E-06, 7.18E-06, 3.13E-06, 4.31E-06, 2.67E-06, 2.12E-06, 5.40E-06, 4.15E-06, 2.51E-06, 1.60E-06, 0.0]

e_high = [4.60E-06, 3.40E-05, 5.25E-05, 4.22E-05, 9.37E-06, 3.13E-05, 1.84E-05, 1.77E-05, 1.46E-05, 1.27E-05, 1.03E-05, 1.63E-05, 7.00E-06, 1.17E-05, 1.32E-05, 1.14E-05, 9.26E-06, 8.24E-06, 1.08E-05, 5.05E-06, 6.98E-06, 2.24E-06, 2.15E-06, 1.20E-06, 2.80E-06, 2.86E-06, 3.41E-06, 3.56E-06, 7.33E-07, 0.0]

lcp = lensingClassFile.lensingClass(0.7, 0.3, 0.7)
F2redshifts = lcp.getRedshiftsFromDLS()
print "Calculating p(z)"
orgList, numList = lcp.organizeRedshifts(F2redshifts)
nlistp = lcp.n_list(numList, orgList, F2redshifts)
#print "nlistp: "
#print nlistp
print "Finished calculating p(z)"
print "numlist:"
print numList
print "orgList:"
print orgList    	


#f = open('nlistp.txt', 'r')
#nlistp = f.read()
#f.close()

#print nlistp

def lnprior(phi):
	omegamat = phi[0]
	omegaDE = phi[1]
	if omegamat < 0.05 or omegamat > 1.0 or omegaDE < 0.0 or omegaDE > 1.0:
		return -inf
	return 0

def lnlike(phi, theta, edat, e_high, e_low, F2red, nlistp):

	#start = time.time()

	omegamat = phi[0]
	omegaDE = phi[1]

	#print("\nshearClacClass:  Omegamat:\t" + str(omegamat) + "  OmegaDE:\t" + str(omegaDE))
	#e_plus = runshear.runshearPS(p[0], p[1]) 
    	shearPSClass = shearPSFile.shearCalcClass(1, 0, 0, 0.673, omegamat, omegaDE, nlistp)
	Pk_l, l = shearPSClass.shearPSPrePrep(F2red)

    	convolutionClassObj = shearPSFile.convolutionClass(1, 0, 0, 0.673, omegamat, omegaDE)
    	e_plus, e_minus, theta = convolutionClassObj.convolutionCalc(Pk_l, l)

	#end = time.time()

	#f = open('stats2.txt', 'w')
	#f.write(str(theta) + '\n')
	#f.close()

	
	sum=0.0
	for i in range(len(edat)):
	    sum += (e_plus[i] - edat[i])**2	 
	   #if e_plus[i] - edat[i] > 0.0:
	#	s_n_2 = e_high[i]**2# + 1.0**2*e_plus[i]**2
	#	print 's_n_2', s_n_2
	#	#sum += pow(edat[i] - e_plus[i], 2.0)/e_high[i]**2 + math.log(2*math.pi*s_n_2)
	#	sum += pow(edat[i] - e_plus[i], 2.0)/s_n_2# + math.log(2*math.pi*s_n_2)

	 #   else:
	#	s_n_2 = e_low[i]**2# + 1.0**2*e_plus[i]**2
	#	print 's_n_2', s_n_2
	#	#sum += pow(edat[i] - e_plus[i], 2.0)/e_low[i]**2 + math.log(2*math.pi*s_n_2)
	#	sum += pow(edat[i] - e_plus[i], 2.0)/s_n_2# + math.log(2*math.pi*s_n_2)

	#return math.exp(-sum)
	return -0.5*sum


def lnprob(p, theta, edat, e_high, e_low, F2red, nlistp):

	lp = lnprior(p)
	if not isfinite(lp):
		return -inf
	return lp + lnlike(p, theta, edat, e_high, e_low, F2red, nlistp)


#import scipy.optimize as opt
#nll = lambda *args: -lnlike(*args)
#result = opt.minimize(nll, [0.3, 0.7],
#                     args=(theta, edat, F2redshifts))
#print result['x']

trueomegamat= 0.315
trueomegaDE = 0.685
trueh = 0.673

##shearPSClass = shearPSFile.shearCalcClass(1, 0, 0, trueh, trueomegamat, trueomegaDE, nlistp)
##Pk_l, l = shearPSClass.shearPSPrePrep(F2redshifts)

##convolutionClassObj = shearPSFile.convolutionClass(1, 0, 0, trueh, trueomegamat, trueomegaDE)
##e_plus, e_minus, theta = convolutionClassObj.convolutionCalc(Pk_l, l)

##print e_plus
##sys.exit()

Nwalkers, Ndim = 8, 2
pos = [[trueomegamat, trueomegaDE] + 0.2*numpy.random.randn(Ndim) for i in range(Nwalkers)]
#pos = [[trueomegamat, trueomegaDE] for i in range(Nwalkers)]
print 'Pos:', pos, '\n\n\n\n\n\n\n\n\n\n\n\n\n\n'
#random.seed(time.time())

#p0 = [0.3, 0.7]# + random.randn(Ndim) for i in range(Nwalker)]

sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnprob, args=(theta, edat, e_high, e_low, F2redshifts, nlistp), threads=16)

pos,prob,state = sampler.run_mcmc(pos, 200)

samples = sampler.chain[:, :, :].reshape((-1, Ndim))

print "pos", pos
print "prob", prob
print "state", state
print "sampler", sampler


res=plt.plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#plt.axhline(alpha_true, color='blue')
plt.axis([0, 25, 0, 1])
plt.savefig("plt.png")
plt.close()

import triangle
#tmp = triangle.corner(sampler.flatchain, labels=['omegamat','omegaDE'], truths=[0.25, 0.75])

fig = triangle.corner(samples)#, labels=["$m$", "$b$"],truths=[0.277, 0.723])
fig.savefig("triangle.png")

fig2 = triangle.corner(samples[:,:], labels=['omegamat', 'omegaDE'], truths = [trueomegamat, trueomegaDE])
fig2.savefig("triangle2.png")

fig3 = triangle.corner(samples[:,:], labels=['omegamat', 'omegaDE'], truths = [trueomegamat, trueomegaDE], extents=[[0.0, 1.0], [0.0, 1.0]])
fig3.savefig("triangle3.png")

#res=plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#axhline(alpha_true, color='blue')

samples[:, 2] = np.exp(samples[:, 2])
omegamatMCMC, omegaDEMCMC = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [33, 50, 68], axis=0)))
print 'omegamatMCMC: ', omegamatMCMC, 'omegaDE: ', omegaDEMCMC

omegamatMCMC, omegaDEMCMC = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
print 'omegamatMCMC: ', omegamatMCMC, 'omegaDE: ', omegaDEMCMC


print "neither here"

sampler.acceptance_fraction

print "nor there"

pymc.Matplot.samplePlot(sampler)

#pymc.gelman_rubin(sampler)

