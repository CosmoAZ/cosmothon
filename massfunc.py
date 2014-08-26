## @package massfunc
# 
# basic mass function calculation
# @note VERY slow!
# [More details.]

import math
import scipy.integrate as integ

import const
import cosmocalcs
import powerspec



## Class that calculates Sheth-Torman mass function
#
class massFunctionST(object):

    ## constructor
    # @param ccalcs      cosmological calculations
    # @param powerspec   power spectrum
    # #param typeLog     if true return dn/dlogM instead of dn/dM
    def __init__(self, ccalcs, powerspec, typeLog=False):
    
        self.ccalcs = ccalcs
        self.powerspec = powerspec
    
        # check ccalcs and power spectrum have same cosmology
        eps = 1e-6
        if ( abs(self.powerspec.omega_m - self.ccalcs.omegamat)>eps ):
            print 'Power spectrum Omega_m =',self.powerspec.omega_m,
            print 'not equal to cosmo calcs Omega_m =',self.ccalcs.omegamat
            return None 
            
        if ( abs(self.powerspec.omega_b - self.ccalcs.omegabaryon)>eps ):
            print 'Power spectrum Omega_b =',self.powerspec.omega_b,
            print 'not equal to cosmo calcs Omega_b =',self.ccalcs.omegabaryon
            return None 
            
        if ( abs(self.powerspec.omega_l - self.ccalcs.omegaDE)>eps ):
            print 'Power spectrum Omega_l =',self.powerspec.omega_l,
            print 'not equal to cosmo calcs Omega_l =',self.ccalcs.omegaDE
            return None 
            
        if ( abs(self.powerspec.h - self.ccalcs.h)>eps ):
            print 'Power spectrum h =',self.powerspec.h,
            print 'not equal to cosmo calcs h =',self.ccalcs.h
            return None


        # set Omega_m
        self.omega_m = self.powerspec.omega_m
        self.h = self.powerspec.h
    
        self.typeLog = typeLog
        
        self.lmstep = 0.05 # step in log of mass
        # step with which to compute derivation of dlnsig/dlnm
        
        # Sheth-Torman mass function parameters
        self.aST = 0.707
        self.AST = 0.3222
        self.QST = 0.3
        
        
    ## Return mass function, either dn/dm or dn/dlogm [in (Mpc/h)^-3 units]
    # @param m    mass in Solar masses/h
    # @param z    redshift
    def getMF(self, m, z=0.):
        
        if (self.typeLog):
            return (self.rhobar0()/m)*self.fST(m,z)*abs(self.dlnsigdlnm(m,z))*math.log(10.)
        else:
            return (self.rhobar0()/(m*m))*self.fST(m,z)*abs(self.dlnsigdlnm(m,z))*math.log(10.)


    ## Return mean density today in units of h^2 M_solar / Mpc^3
    def rhobar0(self):
        return const.rho_crit*self.omega_m
        
        
    ## Differentiatial of ln(variance) by ln(mass) as a function of mass
    # @param m    mass in Solar masses/h
    # @param z    redshift
    def dlnsigdlnm(self, m, z):

        logm2 = math.log(m) + self.lmstep
        m2 = math.exp(logm2)

        sigm1 = math.sqrt(self.sigsqM(m, z))
        sigm2 = math.sqrt(self.sigsqM(m2, z))
        diff = (math.log(sigm2) - math.log(sigm1))/(logm2 - math.log(m))
	
        return diff;

        
    ## Return Sheth-Torman function at mass m and redshift z
    # @param m    mass in Solar masses/h
    # @param z    redshift
    def fST(self, m, z=0.):
    
        # the parameters
        
        # critical density for collapse in spherical collapse model
        deltc = const.DELTA_C

        # variance at R = cuberoot[(3M)/(4PI*RHO)] and redshift z
        sigsq = self.sigsqM(m, z)  
        
        epart = math.exp( -0.5*self.aST*( (deltc*deltc)/sigsq ) )
        endpart = 1. + 1./( self.aST*deltc*deltc/pow(sigsq, self.QST) )

        fst = self.AST*math.sqrt(2.*self.aST/const.pi)*(deltc/math.sqrt(sigsq))*epart*endpart

        return fst
	
	
	## Variance as a function of mass
	# @param m    mass in Solar masses/h
	# @param z    redshift
    def sigsqM(self, m, z):

        # equivalent radius to mass m
        Rcubed = (3*m)/(4.*const.pi*self.rhobar0())
        R = pow(Rcubed, (1./3.));

        # variance of power spectrum at this scale
        sigsq = self.powerspec.getVar(z, R)

        return sigsq
    
    ## Return number of halos in some volume
    # @param fsky    fraction of sky volume covers
    # @param z1      minimum redshift of volume
    # @param z2      maximum redshift of volume
    def numberOfHalos(self, fsky=0.1, z1=0., z2=0.5, mmin=5e13, mmax=1e15):
    
       # integrate volume element
       dVdOmega = self.integrateVolEl(z1, z2)
       vol = 4.*const.pi*fsky*dVdOmega # in Mpc^3
       
       # integrate mass function
       z = (z2 + z1)/2.
       dn = self.integrateMF(mmin, mmax, z) # in (Mpc/h)^3
       
       return dn*(vol*self.h*self.h*self.h)
       
       
    ## Integrate volume element between two redshifts
    # @param zmin
    # @param zmax
    def integrateVolEl(self, zmin, zmax):
        I = integ.quad(self.returndV, zmin, zmax)
        return I[0]
        
        
    ## Return the volume element at redshift in Mpc^3
    # @param z    redshift 
    def returndV(self, z):
        self.ccalcs.setEmissionRedShift(z)
        return self.ccalcs.comovingVolumeElement()
       
       
    ## Integrate mass function between two masses
    # @param mmin    min mass in Solar masses/h
    # @param mmax    max mass in Solar masses/h
    # @param z       redshift of mass function
    def integrateMF(self, mmin, mmax, z=0.):
    
        I = integ.quad(self.getMF, mmin, mmax, args=(z))
        return I[0]
    
