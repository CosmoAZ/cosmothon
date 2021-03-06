import sys
import math
import getopt
import cosmocalcs
import lensingClassFile
import numpy 
import matplotlib.pyplot as plt
from numpy import roll
#import time

def main(argv):
    print ''
    print '=== Running lensing.py program ===='
    """##################    

    lc = lensingClassFile.lensingClass(0)

    lc.calcPowerSpectrum()

    ###################"""

    print "Initializing lc class\n"

    lc = lensingClassFile.lensingClass(0)

    print "Extracting redshifts from DLS data"

    F2redshifts = lc.getRedshiftsFromDLS()

    organizedList, numberedList = lc.organizeRedshifts(F2redshifts)

    print "Printing the organized list after organizeRedshifts function call"
    print "OrgList is: ?", organizedList
    print "Sorting redshift values"

    sortedF2Redshifts = lc.sort_Redshift_Values(organizedList)
    #sortedF2Redshifts = lc.sort_Redshift_Values(F2redshifts[::20])

    Warray = []
    Zarray = []

    for z in sortedF2Redshifts:#[0:len(sortedF2Redshifts) - 1]:

        print "Determining lensing weight function at z: ", z
        
	#start = time.time()
        Wval = lc.lensingWeightFunction(z, sortedF2Redshifts, organizedList, numberedList)
	#end = time.time()

	print "W is: ", Wval

        Warray.append(Wval)
        Zarray.append(z)
    
    ###############################################################
    """
    numberedList = roll(numberedList, 50)

    lc = lensingClassFile.lensingClass(0)

    sortedF2Redshifts = lc.sort_Redshift_Values(organizedList)

    Warray = []
    Zarray = []

    for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:

        print "Determining lensing weight function at z: ", z

        Wval = lc.lensingWeightFunction(z, sortedF2Redshifts, organizedList, numberedList)

        Warray.append(Wval)
        Zarray.append(z)

        f.write(str(z) + '\t' + str(Wval) + '\n')

    plt.plot(Zarray, Warray, 'g--')
    """
    plt.plot(Zarray, Warray, 'g--')
    plt.xlabel('Redshift')
    plt.ylabel('W(z)')
    plt.title('Lensing Weight Function vs. Redshift')
    plt.axis([min(Zarray), max(Zarray), min(Warray), max(Warray)])
    plt.grid(True)
    plt.savefig('WvsZ.png')
    plt.close()

    """plt.plot(Zarray, tempW, 'g--')
    plt.xlabel('Redshift')
    plt.ylabel('\\frac{Da(z\') - Da(z)}{Da(z\')}')
    plt.title('\\frac{Da(z\') - Da(z)}{Da(z\')} vs. Redshift')
    plt.axis([min(Zarray), max(Zarray), min(tempW), max(tempW)])

    plt.grid(True)

    plt.savefig('WunitsvsZ.png')"""


if __name__ == "__main__":
    main(sys.argv[1:])
    

