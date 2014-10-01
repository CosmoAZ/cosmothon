""" Demo of calculating galaxy type fractions as a function of redshift
    
    Input required: galaxy luminosity function read as a table of:
    * first column is redshift
    * following n*3 columns are M*, phi*, alpha for each of the n types
    * each row is a different redshift
    
    Really need to change faint limit according to redshift? and survey mag lim
    
    Need to convert LF absolute magnitude in its observed waveband to the observed apparent magnitude
    in that band given the redshift and galaxy type of the LF
    
    Info needed:
    - waveband of LF
    - transmission of LF
    - spectrum of LF galaxy type
    
    Pre-compute k-corrections for each galaxy type and waveband that LF's are defined for
    Read in those files
    
    LF parameter file will contain all needed info in header like this:
    
    # WAVEBAND OF LF
    # TYPE1, TYPE2, TYPE3
    
    Then code will put together type and waveband and read in correct k-correction

	will need to store the following:
	lumfuncs/
	filters/
	kcorrections/
	spectra/
	

"""

import sys
import getopt
import math
import numpy as np

import galdist # load module that contains Schechter function calculation


def usage(val):
    print 'Usage instructions: '
    print 'typefractions.py -i <lffile> -o <outfile> -m <mfaint=10,mbright=-30>       '
    print ''
    print '<lffile>          file containing luminosity function data                 '
    print '<outfile>         file galaxy type fractions are written to                '
    print '<mfaint,mbright>  magnitude limits to integrate luminosity function between'
    print ''
    print 'Input required: galaxy luminosity function read as a table of:             '
    print '* first column is redshift                                                 '
    print '* following n*3 columns are M*, phi*, alpha for each of the n types        '
    print '* each row is a different redshift                                         '
    print ''
    print 'Empty data values must be entered as nan                                   '
    print ''
    print 'Example:'
    print 'python typefractions.py -i ~/lf.txt -o ~/typefracs.txt          '
    print ''
    sys.exit(val)


def main(argv):
    print ''
    print '=== Running typefractions.py program ===='

    outfile = ''      # file galaxy type fractions are written to
    lffile = ''       # file containing luminosity function data 
    mfaint=10.        # integrate luminosity function to mfaint 
    mbright=-30.      # integrate luminosity function from mbright
    
    try:
        opts, args = getopt.getopt(argv,"ho:i:m:")
    except getopt.GetoptError as err: # if include option that's not there
        usage(2)
      
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt in ("-o"):
            outfile = arg
        elif opt in ("-i"):
            lffile = arg
        elif opt in ("-m"):
            mlims = arg.split(',')
            mfaint = float(mlims[0])
            mbright = float(mlims[1])
            
    print ''
    print 'Galaxy type fractions will be written to file', outfile
    print 'Calculated from luminosity function parameters in', lffile
    print 'Will integrate luminosity function between', mbright,'< M <',mfaint
    print ''
    
    lfdata = np.loadtxt(lffile, comments='#')
    nz = lfdata.shape[0]
    nty = int((lfdata.shape[1] - 1.)/3.)
    print 'Number of redshift luminosity function is defined at =', nz
    print 'Number of galaxy types =', nty
    print ''
  
    integral_array = np.zeros([nz,nty])
    
    for iz in range(nz):
    
        zs = lfdata[iz][0]
        #print 'On redshift', zs,':',
        
        sum_intval = 0.
        for it in range(nty):
            
            icol = 1+it*3
            Mstar = lfdata[iz][icol+0]
            phistar = lfdata[iz][icol+1]
            alpha = lfdata[iz][icol+2]
            
            if (not math.isnan(Mstar)):
                x = galdist.schechterFunction(phistar, Mstar, alpha)
                intval = x.integrate(0., mfaint, mbright)
            else:
                intval = 0.
                
            sum_intval+=intval
            integral_array[iz,it] = intval
            #print intval,
            
        integral_array[iz,:] /= sum_intval
            
        #print ''   
    
    
    for iz in range(nz):
    
        zs = lfdata[iz][0]
        print 'redshift =', zs,'fraction of',
        
        for it in range(nty):
            print 'type',it,'=',integral_array[iz,it],',',
        print ''
    
    # WRITE TO FILE
    #f = open(outfile, 'w')
    
    #for i in range(nk):
        # something
    
    #f.close()
    
    #print '=== File', outfile ,'closed, end of program! ==='
    print ''

if __name__ == "__main__":
    main(sys.argv[1:])
