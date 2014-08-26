## @package powerspec
# 
# linear matter power spectrum with Eisenstein & Hu transfer function and linear LCDM growth function
# @note LCDM only

import math
import scipy.integrate as integ

import const


## Class that calculates Eisenstein & Hu transfer function
#
# @note if Omega_b>Omega_c (particularly if Omega_b>>Omega_c) then every other oscillation has negative 
# transfer function value, therefore only accurate if Omega_b/Omega_c~<2
class transferFunction(object):

    ## constructor
    # @param omega_cdm
    # @param omega_b
    # @param h
    def __init__(self, omega_cdm=0.25, omega_b=0.05, h=1.):

        if (omega_cdm<0. or omega_b<0. or h<0.):
            return -1

        # constructor arguments
        self.omega_c = omega_cdm
        self.omega_b = omega_b
        self.h = h

        # cosmological parameters 
        self.omega_m = self.omega_c + self.omega_b
        H0 = 100.*self.h
        h2 = self.h*self.h
        self.th2p7 = const.T_CMB_K/2.7
        th2p7P4 = self.th2p7*self.th2p7*self.th2p7*self.th2p7


        # Initalization of calculated quantities

        # formula 2 p 606
        self.zeq = 2.50e4*self.omega_m*h2/th2p7P4


        # formula 3 p 607
        # (attention here C=1 : H0 -> H0/C if we use the first formula)
        #  self.keq = sqrt(2.*self.omega_m*H0*H0*self.zeq) / const.cms
        self.keq = 7.46e-2*self.omega_m*h2/(self.th2p7*self.th2p7) # in 1/Mpc units


        # formula 4 p 607
        b1_eq4 = 0.313*pow(self.omega_m*h2, -0.419)*(1. + 0.607*pow(self.omega_m*h2, 0.674))
        b2_eq4 = 0.238*pow(self.omega_m*h2, 0.223)
        self.zd = 1291. * pow(self.omega_m*h2, 0.251) / (1.+0.659*pow(self.omega_m*h2, 0.828)) \
                        * (1. + b1_eq4*pow(self.omega_b*h2, b2_eq4))

                
        # formula 5 page 607    (R = 3*rho_baryon/4*rho_gamma)
        self.Req = 31.5*omega_b*h2/th2p7P4*(1.e3/self.zeq)
        # WARNING: W.Hu's code (tf_fit.c) is in disagreement with the article: zd -> (1+zd)
        self.Rd  = 31.5*omega_b*h2/th2p7P4*(1.e3/self.zd)
        # in tf_fit.c: self.Rd  = 31.5*omega_b*h2 / th2p7P4 * (1.e3/(1.+self.zd));
        

        # formula 6 p 607
        self.s = 2./(3.*self.keq) * math.sqrt(6./self.Req) \
               * math.log( (math.sqrt(1. + self.Rd) + \
                 math.sqrt(self.Rd + self.Req))/(1. + math.sqrt(self.Req)) )


        # formula 7 page 607 (in 1/Mpc units)
        self.ksilk = 1.6*pow(self.omega_b*h2, 0.52)*pow(self.omega_m*h2, 0.73) \
                                                   *(1. + pow(10.4*self.omega_m*h2, -0.95))


        # formulas 11
        a1 = pow(46.9*self.omega_m*h2, 0.670) * (1. + pow(32.1*self.omega_m*h2, -0.532))
        a2 = pow(12.0*self.omega_m*h2, 0.424) * (1. + pow(45.0*self.omega_m*h2, -0.582))
        self.alphac = pow(a1, -self.omega_b/self.omega_m) * pow(a2, -pow(self.omega_b/self.omega_m, 3.))

        # formulas 12
        b1 = 0.944 / (1. + pow(458.*self.omega_m*h2, -0.708))
        b2 = pow(0.395*self.omega_m*h2, -0.0266)
        self.betac = 1. / ( 1. + b1*(pow(self.omega_c/self.omega_m, b2) - 1.) )


        # formula 23 page 610
        self.bnode = 8.41 * pow(self.omega_m*h2, 0.435)


        # formula 15
        # WARNING: W.Hu's code (tf_fit.c) is in disagreement with the article: (1+zeq) -> zeq
        y = (1. + self.zeq)/(1. + self.zd)
        # in tf_fit.c: y = zeq_/(1.+zd_);
        s1py = math.sqrt(1. + y)
        Gy = y*( -6.*s1py + (2. + 3.*y)*math.log((s1py + 1.)/(s1py - 1.)) )
        
        # formula 14 page 608
        self.alphab = 2.07*self.keq*self.s*pow(1. + self.Rd, -3./4.)*Gy

      
        # formula 24 page 610
        self.betab = 0.5 + self.omega_b/self.omega_m   \
                   + (3.-2.*self.omega_b/self.omega_m) * math.sqrt(pow(17.2*self.omega_m*h2,2.) + 1.)

        ### Below here quantites aren't needed for transfer function calculation

        # formula 31 page 612
        self.alphag = 1. \
            - 0.328*math.log(431.*self.omega_m*h2)*self.omega_b/self.omega_m \
            + 0.38*math.log(22.3*self.omega_m*h2)*pow(self.omega_b/self.omega_m, 2.)


        # The approximate value of the sound horizon, formula 26 page 611
        self.sfit = 44.5*math.log(9.83/(self.omega_m*h2)) / math.sqrt(1.+10.*pow(self.omega_b*h2,3./4.))  # Mpc


        # Position of the 1st acoustic peak, formula 25 page 611
        self.kpeak = 5*const.pi/(2.*self.sfit) * (1. + 0.217*self.omega_m*h2)  # 1/Mpc



    ## Return transfer function at wavenumber k
    # @param k   wavenumber in units of h/Mpc
    def returnTF(self, k):

	
        # input k in h/Mpc, but remainder of this function requires it in Mpc^-1
        k = k*self.h 

        # --- CDM
        Tc = self.returnTc(k)

        #--- Baryons
        Tb = self.returnTb(k)

        # --- Total
        T = (self.omega_b/self.omega_m)*Tb + (self.omega_c/self.omega_m)*Tc
        
        return T


    def returnTc(self, k):
    
        # --- CDM
        f = 1./(1. + pow(k*self.s/5.4, 4.) )
        Tc = f*self.T0tilde(k, 1., self.betac) + (1.-f)*self.T0tilde(k, self.alphac, self.betac)

        return Tc


    def returnTb(self, k):

        #--- Baryons
        #  E formula 22 page 610
        ksbnode = k*self.s/self.bnode
        if (ksbnode<0.001): # if ks<bnode so why is this the equality???
            stilde = self.s*ksbnode
        else:
            stilde = self.s / pow(1. + pow(1./ksbnode,3.), 1./3.)

        # spherical Bessel function
        j0kst = 0.
        x = k*stilde
        if (x<0.01):
            j0kst = 1. - x*x/6.*(1. - x*x/20.)
        else:
            j0kst = math.sin(x)/x


        # formula 21 page 610
        Tb = self.T0tilde(k, 1., 1.)/(1. + pow(k*self.s/5.2,2.))
        Tb += self.alphab/(1. + pow(self.betab/(k*self.s), 3.)) * math.exp(-pow(k/self.ksilk, 1.4))
        Tb *= j0kst

        return Tb


    ## Part of fitting function
    # @param k        wavenumber in 1/Mpc units
    # @param alphac   transfer function fit parameter (see eqn 11 in Eisenstein & Hu)
    # @param betac    transfer function fit parameter (see eqn 12 in Eisenstein & Hu)
    def T0tilde(self, k, alphac, betac):

        # WARNING k in NONSTANDARD units of 1/Mpc!
        
        # formula 10 p 608
        q = k/(13.41*self.keq)

        # formula 20 p 610
        C = (14.2/alphac) + 386./(1.+69.9*pow(q, 1.08))

        # formula 19 p 610
        x = math.log(const.e + 1.8*betac*q)
        return x / (x + C*q*q)


## Growth function in Lambda-CDM universe (w=-1)
# 
class growthFunctionLCDM(object):

    ## constructor
    # @param omega_m
    # @param omega_l
    def __init__(self, omega_m=0.25, omega_l=0.75):
    
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_k = 1. - self.omega_m - self.omega_l
        self.D0 = 1. # temporarily set to 1
        self.D0 = self.getGrowth(0.)
        
        
    ## Return growth at redshift z (normalized to 1 at z=0)
    # @param z    redshift
    def getGrowth(self, z):
        
        omz = self.omegaMz(z) # Omega_m(z)
        olz = self.omegaLz(z) # Omega_l(z)
    
        D = (2.5*omz)/( (1.+z)*( pow(omz,4./7.) - olz + (1. + 0.5*omz)*(1. + olz/70.)) )
        return D/self.D0
    
    
    # Return Omega_m at redshift z
    # @param z    redshift
    def omegaMz(self, z):
        
        z2 = (1. + z)*(1. + z) # just (1+z)^2
        omz = (self.omega_m*z2*(1. + z))/(self.omega_l + self.omega_k*z2 + self.omega_m*z2*(1. + z))
        return omz
    
    
    # Return Omega_l at redshift z
    # @param z    redshift
    def omegaLz(self, z):
    
        z2 = (1. + z)*(1. + z) # just (1+z)^2
        olz = self.omega_l/(self.omega_l + self.omega_k*z2 + self.omega_m*z2*(1. + z))
        return olz
        
        
## Power spectrum calculated using LCDM growth function and Eisenstein & Hu transfer function
# 
#  Normalised with the amplitude of primordial curvature fluctuations at the pivot scale k =0.05 Mpc^-1
class powerSpectrum(object):

    ## Constructor
    # @param transfunc    transfer function
    # @param grofunc      growth function
    # @param ns           spectral index
    # @param DelRsq       amplitude of primordial curvature fluctuations at the pivot scale
    def __init__(self, transfunc, grofunc, ns=1., DelRsq=2.1e-9):

        self.transfunc = transfunc
        self.grofunc = grofunc
        self.ns = ns
        self.DelRsq = DelRsq
        
        if (self.transfunc.omega_m != self.grofunc.omega_m):
            print 'Transfer function Omega_m =',self.transfunc.omega_m,
            print 'not equal to growth function Omega_m =',self.grofunc.omega_m
            return None
         
        print 'Transfer function has parameters Omega_m =',self.transfunc.omega_m,'and h =',self.transfunc.h
        print 'Growth function has parameters Omega_m =',self.grofunc.omega_m,
        print 'and Omega_L =',self.grofunc.omega_l
         
        self.omega_m = self.transfunc.omega_m
        self.h = self.transfunc.h
        self.omega_b = self.transfunc.omega_b
        self.omega_l=self.grofunc.omega_l
        
        self.kpivot = 0.05
        self.Apiv = 8.0967605e7
        
           
    ## Return value of power spectrum at wavenumber k and redshift z
    # @param k    wavenumber k in h/Mpc
    # @param z    redshift
    def getPk(self, k, z=0.):
    
        tf = (self.transfunc).returnTF(k)
        power = self.power()
        knorm = self.kNorm(k)
        h4 = self.hpower4()
        d1 = (self.grofunc).getGrowth(z)
        powfactor = self.powFactor()
        
        # dimensionless power spectrum Delta_k^2
        DeltaSqK = self.Apiv*pow(knorm, power)*tf*tf*powfactor*self.DelRsq*d1*d1/h4
        
        # P(k)
        pspec = (2.*const.pi*const.pi)*DeltaSqK/(k*k*k)
        
        return pspec

    
    # Return spectral index plus 3
    def power(self):
        return self.ns + 3.
        
    
    # Return k normalized to pivot scale
    def kNorm(self, k):
        return (self.h*k)/self.kpivot
        
    
    # Return Hubble constant to the fourth power
    def hpower4(self):
        return self.h*self.h*self.h*self.h
        
        
    # Return growth at redshift zero divided by Omega_m all squared \f$ (D(z=0)/\Omega_m)^2 \f$
    def powFactor(self):
        d10 = (self.grofunc).D0
        return (d10/self.omega_m)*(d10/self.omega_m)
        
        
    ## Return power spectrum variance \f$ \sigma^2(R)=1/2\pi \int dk/k k^3 P(k;z) W^2(kR) \f$
    # @param z       redshift of power spectrum
    # @param R       scale in Mpc/h variance measured over (top hat filter size)
    # @param kmin    integrate power spectrum integrand from this wavenumber
    # @param kmax    integrate power spectrum integrand to this wavenumber
    def getVar(self, z=0., R=8., kmin=1e-6, kmax=1000.):
    
        I = integ.quad(self.getIntegrand, kmin, kmax, args=(z,R))
        return I[0]
    
    
    ## Return integrand for power spectrum variance calculation
    # @param k       wavenumber k in h/Mpc
    # @param z       redshift of power spectrum
    # @param R       scale in Mpc/h variance measured over (top hat filter size)
    def getIntegrand(self, k, z, R):
        return k*k*self.getPk(k, z)*self.Wk(k*R)*self.Wk(k*R)/(2.*const.pi*const.pi)
    
    
    ## Return top hat filter
    # @param kR    wavenumber k in h/Mpc multiplied by size of top hat filter in Mpc/h
    def Wk(self, kR):
        if (kR<0.1):
            return -0.1*kR*kR + 1.
        else:
            return 3.*(math.sin(kR) - kR*math.cos(kR))/(kR*kR*kR)
    
