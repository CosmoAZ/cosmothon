""" Demo of matter power spectrum calculation
    
    The matter power spectrum has 3 components:
    * Initial power spectrum P_i(k) = Ak^n
    * Transfer function: T(k)
    * Growth function: D(z)
    
    Matter power spectrum at redshift z, wavenumber k (in h/Mpc units):
    P_m(k,z) = D(z)^2*P_i(k)*T(k)^2
    
"""

import sys
import getopt
import math
import powerspec # load module that contains power spectrum calcs


def usage(val):
    print 'Usage instructions: '
    print 'matterpowerspec.py -o <outfile> -z <redshift> -k <kmin,kmax> -c <om,ol,ob,h>      '
    print ''
    print '<outfile>     file matter power spectrum written to                               '
    print '<redshift>    redshift to matter power spectrum at                                '
    print '<kmin,kmax>   min and max wavenumbers to calculate matter power spectrum between  '
    print '<om,ol,ob,h>  cosmological parameters to use in calculation:                      '
    print '              om=Omega_m: total matter density (dark matter+baryonic)             '
    print '              ol=Omega_lambda: dark energy density                                '
    print '              ob=Omega_b: baryonic matter density                                 '
    print '              h=unitless Hubble constant                                          '
    print ''
    print 'Example:'
    print 'python matterpowerspec.py -o ~/ps.txt -z 1. -k 0.0001,1 -c 0.3,0.7,0.05,0.72      '
    print ''
    sys.exit(val)


def main(argv):
    print ''
    print '=== Running matterpowerspec.py program ===='

    outfile = ''      # file matter power spectrum written to
    redshift = 0.     # redshift of matter power spectrum
    kmin = 0.001      # min wavenumber of matter power spectrum 
    kmax = 1.         # max wavenumber of matter power spectrum
    omega_m = 0.3     # total matter density (dark matter+baryonic)  
    omega_l = 0.7     # dark energy density 
    omega_b = 0.05    # baryonic matter density 
    h = 0.72          # unitless Hubble constant  
    
    try:
        opts, args = getopt.getopt(argv,"ho:z:k:c:")
    except getopt.GetoptError as err: # if include option that's not there
        usage(2)
      
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt in ("-o"):
            outfile = arg
        elif opt in ("-z"):
            redshift = float(arg)
        elif opt in ("-k"):
            kvals = arg.split(',')
            kmin = float(kvals[0])
            kmax = float(kvals[1])
        elif opt in ("-c"):
            pars = arg.split(',')
            omega_m = float(pars[0])
            omega_l = float(pars[1])
            omega_b = float(pars[2])
            h = float(pars[3])
            
    print ''
    print 'Matter power spectrum will be written to file', outfile
    print 'It will be calculated at redshift', redshift ,'and wavenumbers', kmin ,'h/Mpc to', kmax ,'h/Mpc' 
    print 'For the following cosmological parameters:'
    print 'Omega_m =', omega_m
    print 'Omega_l =', omega_l
    print 'Omega_b =', omega_b
    print 'h =', h
    print '(Dark energy equation of state w0=-1, wa=0)'
    print ''


    # GROWTH FUNCTION
    gro = powerspec.growthFunctionLCDM(omega_m, omega_l)

    
    # TRANSFER FUNCTION
    omega_c = omega_m - omega_b  # cold dark matter density
    tf = powerspec.transferFunction(omega_c, omega_b, h)


    # TOTAL POWER SPECTRUM
    
    ns = 0.96  # spectral index
    DelRsq = 2.1e-9  # amplitude of primordial curvature fluctuations at the pivot scale
    # Initial power spectrum P_i(k) ~ DelRsq*k^ns
    
    ps = powerspec.powerSpectrum(tf, gro, ns, DelRsq)
    print 'Amplitude of matter fluctuations on 8Mpc/h scale: sigma_8 =',math.sqrt(ps.getVar())
    print ''

    # WRITE TO FILE
    f = open(outfile, 'w')
    
    # calculate power spectrum on log scale
    nk = 1000 # number of points to calculate power spectrum at between kmin, kmax
    lkmin = math.log(kmin)
    lkmax = math.log(kmax)
    dlogk = (lkmax - lkmin)/(nk - 1.)

    for i in range(nk):
    
        k = math.exp(lkmin + i*dlogk) # in h/Mpc units
        Pk = ps.getPk(k) # in (Mpc/h)^3 units
        #wth = ps.Wk(k*8)
    
        f.write(str(k) + "  " + str(Pk) + "\n")
    
    f.close()
    
    print '=== File', outfile ,'closed, end of program! ==='
    print ''

if __name__ == "__main__":
    main(sys.argv[1:])
