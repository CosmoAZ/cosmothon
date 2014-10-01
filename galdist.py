## @package galdist
# 
# this module hold calculations for the galaxy distribution
#
# [More details.]

import math
import os.path
import numpy as np
import scipy.integrate as integ


import const
import cosmocalcs


## Class for computing n(z)
#
# ## default n(z) z0 parameter matches DLS median redshift (=0.78)
#
# Last modified: 30 July 2014 by AA
class nOfz(object):

    ## The constructor
    # @param z0            analytic n(z) function parameter (related to median z)
    # @param alpha         analytic n(z) function parameter
    # @param beta          analytic n(z) function paramete
    def __init__(self, z0=0.41, alpha=2., beta=1.2): 
        self.z0 = z0
        self.alpha = alpha
        self.beta = beta


    ## Return value of n(z) function
    # @param z    redshift
    def getNofZ(self, z):
        return pow(z, self.alpha)*math.exp(-pow(z/self.z0, self.beta))
    
    
    ## Return redshift at maximum n(z)
    def zAtMaxNofz(self):
        return  pow((self.alpha*pow(self.z0,self.beta))/self.beta,1./self.beta)
        
    
    ## Integrate the n(z) between z1 and z2
    # @param z1
    # @param z2
    def integrate(self, z1, z2):
        
        I = integ.quad(self.getNofZ, z1, z2)
        return I[0]
    
    
## Function to return schechter functional form
#
# Last modified: 30 Sept 2014 by AA
def schechter(mag, phistar, mstar, alpha):
    y = pow(10.,-0.4*(mag - mstar))
    return 0.4*math.log(10.)*phistar*pow(y,1+alpha)*math.exp(-y)


## Class for computing schechter function
#
# Last modified: 8 July 2014 by AA
class schechterFunction(object):

    ## The constructor
    # @param phistar    function parameter (units h^3/Mpc^3)
    # @param mstar      function parameter (magnitude units)
    # @param alpha      function parameter (no units)
    def __init__(self, phistar=1.49e-2, mstar=-20.44, alpha=-1.05, mlimit=24.):
    
        self.phistar = phistar
        self.mstar = mstar
        self.alpha = alpha
        self.mlimit = mlimit
        self.ccalc = cosmocalcs.cosmologyCalculator()
        
        
    ## Return value of schechter function at absolute magnitude M and redshift z
    # @param mag    absolute magnitude
    # @param z      redshift (if want to evoluve mstar linearly)
    def getphi(self, mag, z=0.):
        
        magz = mag - z # this controls evolution in mstar
        y = pow(10.,-0.4*(magz - self.mstar))
        return schechter(magz, self.phistar, self.mstar, self.alpha)


    ## Return approximate number density at redshift given the magnitude limit
    # @param z       redshift
    # @param mlim    magnitdue limit
    def getnumberDensity(self, z, mlim):

        # approx conversion from magnitude limit to equivalent absolute magnitude at z
        self.ccalc.setEmissionRedShift(z)
        dL = self.ccalc.LuminosityDistanceMpc()
        Mlim = mlim -  5.*math.log10(dL) - 25.  # faintest ABSOLUTE mag observable at z
        #print 'MLIM',Mlim,mlim,z, dL
        return self.integrate(z, Mlim)
        

    ## Integrate the schechter function at absolute magnitude M and redshift z
    # @param z        redshift
    # @param mmax     faintest absolute magnitude
    # @param mfaint   brightest absolute magnitude
    def integrate(self, z=0., mfaint=10, mbright=-30.):
        
        I = integ.quad(self.getphi, mbright, mfaint, args=(z))
        return I[0]

