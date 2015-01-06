import sys
import shearPSFile

def main(argv):
    print ''
    print '=== Running nonlinearPS.py program ====\n'
    print argv

    i = int(argv[3])
    j = int(argv[5])
    shearBins = int(argv[1])
    h = float(argv[7])
    omegamat = float(argv[9])
    omegaDE = float(argv[11])   

    """for opt, arg in opts:
        if opt == '-i':
            i = arg
        elif opt == '-j':
            j = arg
        elif opt in ("-n"):
            shearBins = arg
        elif opt in ("-h"):
            print "Run program with: python nonlinearPS.py -n <number of shear bins> -i <first bin to calculate cross-spectra spectrum> -j <second bin to calculate cross-spectra spectrum>"

    if i == '' or j == '' or shearBins == '':
	print "i, j, or shearBins not found\n"""

    print "Number of bins is: ", shearBins, " i is: ", i, " j is: ", j

    shearPSClass = shearPSFile.shearCalcClass(shearBins, i, j, h, omegamat, omegaDE)
    Pk_l = shearPSClass.shearPSCalc()
    #shearPSClass.printPSResults()

    #Convolution calculation
    convolutionClassObj = shearPSFile.convolutionClass(shearBins, i, j, h, omegamat, omegaDE)
    e_plus, e_minus, theta = convolutionClassObj.convolutionCalc(Pk, l)
    convolutionClassObj.printCCresults(e_plus, e_minus, theta)



if __name__ == "__main__":
    main(sys.argv[1:])
    
