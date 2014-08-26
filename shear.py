# calculate shear due to circularly symmetric profile
# need to rewrite these to work with new cosmology class

import math
import cmath
import const
import cosmocalcs
import numpy as np

def sisshear(theta, sigma, zcluster, zsource, cc):
    """
    return the tangential shear of a source galaxy at z=zsource expected 
    assuming a SIS profile with velocity dispersion sigma when the galaxy is at 
    radial angular distance=theta (in arcmin) from the center of a galaxy
    cluster at z=zcluster
    
    (number, number, number, number, cosmologyCalculator) -> number
    (list|array|number, number, number, list|array|number, cosmologyCalculator) 
                                                          -> list|array|number
    
    theta:       angular distance in arcmin
    sigma:       velocity dispersion of SIS model in km/s
    zcluster:    redshift of cluster
    zsource:     redshift of background source
    cc:          cosmology calculator class
    
    Assumes a circularly symmetric SIS profile, and a cosmological parameters
    are set in cosmology calculator class

    Assumes flat universe and LCDM
    
    References to Schneider refer to Saas Fee notes:
    http://www.astro.uni-bonn.de/~peter/SaasFee.html
   
    """

    arcmin_per_radian = 60.*180./const.pi
    
    # unpack cosmological model, only implemented OmegaM for now
    if (cc.isClosed() or cc.isOpen()):
        print 'SIS shear model not implemented for non-flat universes yet!',
        print 'OmegaK = ',cc.omegaCurv
        return -99
    
    if (not cc.isLambda):
        print 'SIS shear model not implemented for non-LCDM universes yet!',
        print 'w0 = ',cc.wX,'wa =',cc.wXa
        return -99
    
    OmegaM = cc.omegamat
    OmegaQ = 1. - OmegaM
    w0 = -1.
    wa = 0.
    h = 1.
    cc.sethubble(h) # all in h units
    
    # return Hubble distance c/H0 in Mpc/h
    Dh = cc.DH()
    
    # LENS QUANITITES
    cc.setEmissionRedShift(zcluster)
    # transverse comoving distance to lens Mpc/h
    Dm_lens = cc.TransComovDistanceMpc() 
    
    if (isinstance(zsource, (int, long, float, complex))): # if single z
    
        # SOURCE QUANTITIES
        cc.setEmissionRedShift(zsource)
        # angular diameter distance etc to source Mpc/h
        Ds = cc.AngularDiameterDistanceMpc()
        # transverse comoving distance to source Mpc/h
        Dm_source = cc.TransComovDistanceMpc() 
    
        # Angular diameter distance between lens and source
        sqrt1 = Dm_source*math.sqrt(1.+cc.omegaCurv*Dm_lens*Dm_lens/(Dh*Dh))
        sqrt2 = Dm_lens*math.sqrt(1.+cc.omegaCurv*Dm_source*Dm_source/(Dh*Dh))
        Dds = (sqrt1 - sqrt2)/(1.+zsource)
        
    else: # if list/array of z
        Ds = np.zeros([len(zsource)])
        Dds = np.zeros([len(zsource)])
        theta_hold = np.zeros([len(zsource)])
        i = 0
        for zs in zsource:
        
            # SOURCE QUANTITIES
            cc.setEmissionRedShift(zs)
            # angular diameter distance etc to source Mpc/h
            Ds[i] = cc.AngularDiameterDistanceMpc()
            # transverse comoving distance to source Mpc/h
            Dm_source = cc.TransComovDistanceMpc() 
    
            # Angular diameter distance between lens and source
            sqrt1 = Dm_source*math.sqrt(1.+cc.omegaCurv*Dm_lens*Dm_lens/(Dh*Dh))
            sqrt2 = Dm_lens*math.sqrt(1.+cc.omegaCurv*Dm_source*Dm_source/(Dh*Dh))
            Dds[i] = (sqrt1 - sqrt2)/(1.+zs)
            
            theta_hold[i] = theta[i]
            
            i+=1
        theta = theta_hold
    
    # Schneider eqn 52: Einstein radius
    theta_E = 4.*const.pi*(sigma/const.cKms)*(sigma/const.cKms)*Dds/Ds # radians
    theta_E *= arcmin_per_radian # arminutes
    
    # Schneider eqn 53: 
    kappa = theta_E/(2.*theta)
    mod_gamma = theta_E/(2.*theta)
    
    return mod_gamma, kappa
    

def nfwshear(theta, nfwpars, zcluster, zsource, cc):
    """
    return the tangential shear of a source galaxy at z=zsource expected 
    assuming a NFW profile with parameters=nfwpars when the galaxy is at 
    radial angular distance=theta (in arcmin) from the center of a galaxy
    cluster at z=zcluster
    
    Assumes a circularly symmetric NFW profile with parameters in dictionary
    variable nfwpars, and a cosmology 
    
    (number, dict, number, number, cosmologyCalculator) -> number
    (list|array|number, dict, number, list|array|number, cosmologyCalculator) 
                                                        -> list|array|number
    
    theta:       angular distance in arcmin
    nfwpars:     dictionary containing M200 and concentration
    zcluster:    redshift of cluster
    zsource:     redshift of background source
    cc:          cosmology calculator class

    Assumes flat universe and LCDM
    
    Checked against Sarah's matlab code
    
    
    """
    
    # unpack cosmological model, only implemented OmegaM for now
    if (cc.isClosed() or cc.isOpen()):
        print 'SIS shear model not implemented for non-flat universes yet!',
        print 'OmegaK = ',cc.omegaCurv
        return -99
    
    if (not cc.isLambda):
        print 'SIS shear model not implemented for non-LCDM universes yet!',
        print 'w0 = ',cc.wX,'wa =',cc.wXa
        return -99
    
    
    # unpack NFW model
    m200 = nfwpars["M200"] 
    conc = nfwpars["c"]
    #print 'NFW model: m200 =',m200,' c =',conc

    # unpack cosmological model, only implemented OmegaM for now
    OmegaM = cc.omegamat
    OmegaQ = 1. - OmegaM
    w0 = -1.
    wa = 0.
    h = 1.
    cc.sethubble(h) # all in h units
    
    # mean density at cluster redshift
    z_plus1_cubed = (1.+zcluster)*(1.+zcluster)*(1.+zcluster)
    rho_m_zclust = const.rho_crit*OmegaM*z_plus1_cubed
    #print 'rho_crit =',const.rho_crit,', zcluster =',zcluster
    #print 'rho_m(z_clust) =',rho_m_zclust
    
    # return Hubble distance c/H0 in Mpc/h
    Dh = cc.DH()
    
    # LENS QUANITITES
    cc.setEmissionRedShift(zcluster)
    # angular diameter distance etc to lens Mpc/h 
    Dd = cc.AngularDiameterDistanceMpc()
    # transverse comoving distance to lens Mpc/h
    Dm_lens = cc.TransComovDistanceMpc() 
    
    # NFW model
    delta_c = 200./3.*conc*conc*conc / (math.log(1.+conc) - conc/(1.+conc))
    rho_s = delta_c*rho_m_zclust
    r200 = math.pow(m200/(200.*rho_m_zclust*4./3.*const.pi),1./3.)
    theta200 = (r200/Dd) *180./const.pi*60 # r200 in theta in arcmin
    r_s = r200/conc
    #print 'NFW model parameters: delta_c =',delta_c,' rho_s =',rho_s,'r200 =',r200,
    #print 'theta200 =',theta200,'r_s =',r_s
    
    wasNum = False
    if (isinstance(zsource, (int, long, float, complex))): # if single z
    
        # SOURCE QUANTITIES
        cc.setEmissionRedShift(zsource)
        # angular diameter distance etc to source Mpc/h (don't actually need this)
        Ds = cc.AngularDiameterDistanceMpc()
        # transverse comoving distance to source Mpc/h
        Dm_source = cc.TransComovDistanceMpc() 
    
    
        # Angular diameter distance between lens and source
        sqrt1 = Dm_source*math.sqrt(1.+cc.omegaCurv*Dm_lens*Dm_lens/(Dh*Dh))
        sqrt2 = Dm_lens*math.sqrt(1.+cc.omegaCurv*Dm_source*Dm_source/(Dh*Dh))
        Dds = (sqrt1 - sqrt2)/(1.+zsource)
        
        wasNum = True

    else: # if list/array of z
    
        Ds = np.zeros([len(zsource)])
        Dds = np.zeros([len(zsource)])
        theta_hold = np.zeros([len(zsource)])
        i = 0
        for zs in zsource:
        
            # SOURCE QUANTITIES
            cc.setEmissionRedShift(zs)
            # angular diameter distance etc to source Mpc/h
            Ds[i] = cc.AngularDiameterDistanceMpc()
            # transverse comoving distance to source Mpc/h
            Dm_source = cc.TransComovDistanceMpc() 
    
            # Angular diameter distance between lens and source
            sqrt1 = Dm_source*math.sqrt(1.+cc.omegaCurv*Dm_lens*Dm_lens/(Dh*Dh))
            sqrt2 = Dm_lens*math.sqrt(1.+cc.omegaCurv*Dm_source*Dm_source/(Dh*Dh))
            Dds[i] = (sqrt1 - sqrt2)/(1.+zs)
            
            theta_hold[i] = theta[i]
            
            i+=1
        theta = theta_hold

    
    # sigma crit in Msolar/Mpc^2
    sig_crit = cosmocalcs.Sigma_crit_Msolar_Mpcsq(zcluster, zsource, cc)
    #sig_crit = cosmocalcs.Sigma_crit(zcluster, zsource, 100., OmegaM, OmegaQ, w0, wa)/const.msolarKg
    #print 'Sigma_crit =',sig_crit

    # dimensionless lens parameter
    x = theta*(1./60.)*(const.pi/180.)*(Dd/r_s)
    #print 'x =',x
    
    #print 'If x=1: kappa =',2.*r_s*rho_s/sig_crit/3.
    #print 'If x=1: mod_gamma =', r_s*rho_s/sig_crit*(10./3. + 4.*math.log(0.5))
    
    
    if (wasNum): 
        xvals = np.zeros([1])
        kappa = np.zeros([1])
        mod_gamma = np.zeros([1])
        xvals[0] = x
        sc = sig_crit
        sig_crit = np.zeros([1])
        sig_crit[0] = sc
    else:
        xvals = x
        kappa = np.zeros([len(x)])
        mod_gamma = np.zeros([len(x)])
        #print type(xvals)
        
    i=0
    for x in xvals:

        # kappa and gamma
        if (x<1.):
            #print 'x<1'
            kappa[i] =  2.*r_s*rho_s/sig_crit[i]/(x*x-1.)*(1. - 2./math.sqrt(1 - x*x) \
                       *math.atanh(math.sqrt( (1-x)/(1+x) )) )
            mod_gamma[i] = r_s*rho_s/sig_crit[i]*(4./(x*x)*math.log(x/2.) - 2./(x*x-1) \
                       + 4.*math.atanh(math.sqrt((1-x)/(1+x)))*(2 - 3*x*x) \
                       / (x*x*math.sqrt((1-x*x)*(1-x*x)*(1-x*x))) )
        elif (x == 1.):
            #print 'x=1'
            kappa[i] = 2.*r_s*rho_s/sig_crit[i]/3.
            mod_gamma[i] = r_s*rho_s/sig_crit[i]*(10./3. + 4.*math.log(0.5))
        elif (x>1.):
            #print 'x>1'
            k = 2.*r_s*rho_s/sig_crit[i]/(x*x - 1.)*(1. - 2./math.sqrt(x*x - 1) \
                      *math.atan(math.sqrt((x-1.)/(1.+x) )) )
            kappa[i] = k
            mod_gamma_atan = (r_s*rho_s/sig_crit[i]*(4./(x*x)*math.log(x/2.) - 2./(x*x-1) \
                            + 4.*math.atan(math.sqrt((x-1.)/(1.+x)))*(2. - 3.*x*x)  \
                            / (x*x*math.sqrt(pow(x*x-1.,3))) ) )
        
            # something doesn't seem to work right with atan, so use atanh
            # and take real part (imaginary parts of the below are zero anyway).
            mod_gamma_complex = r_s*rho_s/sig_crit[i]*(4./(x*x)*math.log(x/2.) - 2./(x*x-1)  \
                               + 4.*cmath.atanh(cmath.sqrt((1.-x)/(1.+x)))*(2. - 3.*x*x) \
                               / ( x*x*cmath.exp(1.5*cmath.log(1.-x*x)))  )
        
            if (abs(mod_gamma_complex.imag))>1e-6:
                print 'WARNING! Complex part might be large',abs(mod_gamma_complex.imag)
            mod_gamma[i] = mod_gamma_complex.real
        i+=1         
        #print mod_gamma,type(mod_gamma),mod_gamma_atan,type(mod_gamma_atan)
         
    if (wasNum):
        mod_gamma = mod_gamma[0]
        kappa = kappa[0]
        
        # NB. this agrees with numerical integration of kappa values too
    return mod_gamma, kappa
