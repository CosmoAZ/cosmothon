import sys
import math
import getopt
import cosmocalcs
import lensingClassFile

def main(argv):
    print ''
    print '=== Running lensing.py program ===='

    z = 1


    print "Initializing lc class\n"

    lc = lensingClassFile.lensingClass(1.)


    print "Extracting redshifts from DLS data"

    F2redshifts = lc.getRedshiftsFromDLS()
    

    sortedF2Redshifts = lc.sort_Redshift_Values(F2redshifts)

    for z in sortedF2Redshifts:

        print "Determining lensing weight function\n"

        wat = lc.lensingWeightFunction(z)

        print wat

        


if __name__ == "__main__":
    main(sys.argv[1:])
    

