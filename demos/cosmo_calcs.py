""" Demo of cosmology calculations 

    How to use:
    cosmocalcs.cosmologyCalculator()

    Demonstrates calculation of angular diameter distance and luminosity distance
    by reproducing data for Figures 2 and 3 in Hogg 1999
    
    Also demonstrates calculating angular diameter distance and luminosity distance
    for various dark energy equation of state parameters
"""

import sys
import getopt
import math
import cosmocalcs # load module that contains cosmology calcs


def usage(val):
    print 'Usage instructions: '
    print 'cosmo_calcs.py -o <outroot>                        '
    print ''
    print '<outfile>     root name of file outputs go to      '
    print ''
    print 'Example:'
    print 'python cosmo_calcs.py -o ~/cosmos                  '
    print ''
    sys.exit(val)


def main(argv):
    print ''
    print '=== Running cosmo_calcs.py program ===='
    
    outroot = ''      # root name of file outputs go to 

    try:
        opts, args = getopt.getopt(argv,"ho:")
    except getopt.GetoptError as err: # if include option that's not there
        usage(2)
      
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt in ("-o"):
            outroot = arg

            
    print ''
    print 'Output will be written to file', outroot
    print ''


    # set up default parameters
    h=0.7  # hubble parameter
    wX=-1. # dark energy equation of state today
    wXa=0. # dark energy equation of state evolution


    # Omega_m and Omega_L to match examples in Hogg 1999
    omegams = [1, 0.05, 0.2] 
    omegals = [0., 0., 0.8]



    # Calculate angular diameter distance and luminosity distance to match
    # plots in Hogg 1999

    cosmoEdS = cosmocalcs.cosmologyCalculator(h,omegams[0],omegals[0],wX,wXa) # Einstein de Sitter
    cosmoLowDens = cosmocalcs.cosmologyCalculator(h,omegams[1],omegals[1],wX,wXa)
    cosmoHighLam = cosmocalcs.cosmologyCalculator(h,omegams[2],omegals[2],wX,wXa)
    
    
    # set radiation component of universe (i.e. relativistic constituents like photons, neutrinos) to zero
    # (really size of this component is like 10^-5 so negligible)
    cosmoEdS.setOmegaRadiation(0.)
    cosmoLowDens.setOmegaRadiation(0.)
    cosmoHighLam.setOmegaRadiation(0.)

    
    # redshift range for calculation
    zmin = 0.
    zmax = 5.
    npt = 100
    dz = (zmax-zmin)/(npt-1)

    
    # Open file to write to 
    outfile = outroot + '_hogg.txt'
    f = open(outfile,'w')

    for i in range(npt):

        z = zmin + i*dz
    
        # Set emission redshift means the calculation can be made incrementally from the 
        # previous redshift that was set. This makes calculation faster IF it is repeated for 
        # small steps increasing in redshift
        cosmoEdS.setEmissionRedShift(z)
        cosmoLowDens.setEmissionRedShift(z)
        cosmoHighLam.setEmissionRedShift(z)
        
        # note that this is returning the UNITLESS distances:
        # In Hogg 1999 notation, D_A/D_H and D_L/D_H are being returned
        # (D_H = c/H_0)
        
        # To return these distances in Mpc units use:
        # AngularDiameterDistanceMpc()
        # LuminosityDistanceMpc()
        
        da1 = cosmoEdS.AngularDiameterDistance()
        dl1 = cosmoEdS.LuminosityDistance()
    
        da2 = cosmoLowDens.AngularDiameterDistance()
        dl2 = cosmoLowDens.LuminosityDistance()
    
        da3 = cosmoHighLam.AngularDiameterDistance()
        dl3 = cosmoHighLam.LuminosityDistance()

        # write to file
        row =  str(z) + "  " + str(da1) + "  " + str(dl1) + "  " 
        row += str(da2) + "  " + str(dl2) + "  "
        row += str(da3) + "  " + str(dl3) + "\n"
        f.write(row)
    
    f.close()


    # Take 'high lambda' parameters but vary w0 parameter, calculate angular 
    # diameter distance and luminosity distance 
    now = 3
    wmin=-1.1
    wmax=-0.9
    dw=(wmax-wmin)/(now-1)

    for iw in range(now):

        wv = wmin + iw*dw
        cc = cosmocalcs.cosmologyCalculator(h,omegams[2],omegals[2])
        cc.setOmegaRadiation(0.)
        cc.setDarkEnergyEoS(wv) # this only sets constant w (no evolution)
        #cc.print_pars()
    
        outfile = outroot + "_w=" + str(wv) + ".txt"
        f = open(outfile,'w')

        for i in range(npt):

            z = zmin + i*dz
    
            cc.setEmissionRedShift(z)

            da = cc.AngularDiameterDistance()
            dl = cc.LuminosityDistance()

            row =  str(z) + "  " + str(da) + "  " + str(dl) + "\n" #+ "  " + str(tmp) + "\n"
            f.write(row)
    
        f.close()
    

    # Take 'high lambda' parameters but vary wa parameter, calculate angular 
    # diameter distance and luminosity distance 
    wamin=-0.1
    wamax=0.1
    dwa=(wamax-wamin)/(now-1);

    for iw in range(now):

        wav = wamin + iw*dwa;
        cc = cosmocalcs.cosmologyCalculator(h,omegams[2],omegals[2])
        cc.setOmegaRadiation(0.)
        cc.setDarkEnergyEoS(wX, wav) # two arguments sets w with evolution
    
        outfile = outroot + "_wa=" + str(wav) + ".txt"
        f = open(outfile,'w')

        for i in range(npt):

            z = zmin + i*dz
    
            cc.setEmissionRedShift(z)

            da = cc.AngularDiameterDistance()
            dl = cc.LuminosityDistance()

            row =  str(z) + "  " + str(da) + "  " + str(dl) + "\n" #+ "  " + str(tmp) + "\n"
            f.write(row)
    
        f.close()
    
    
    print '=== end of program! ==='
    print ''

if __name__ == "__main__":
    main(sys.argv[1:])
