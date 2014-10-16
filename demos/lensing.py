import sys
import math
import getopt
import cosmocalcs
import lensingClassFile

def main(argv):
    print ''
    print '=== Running lensing.py program ===='

    z = 1

    f = open('results.txt', 'w')

    print "Initializing lc class\n"

    lc = lensingClassFile.lensingClass(1.)


    print "Extracting redshifts from DLS data"

    F2redshifts = lc.getRedshiftsFromDLS()

    print "Sorting redshift values"

    sortedF2Redshifts = lc.sort_Redshift_Values(F2redshifts)

    i = 1

    for z in sortedF2Redshifts:

        print "Determining lensing weight function at z: ", z
        
        wat = lc.lensingWeightFunction(z, sortedF2Redshifts)
        
        f.write(sortedF2Redshifts[i] + '\t' + wat + '\n')

        i += 1
        


if __name__ == "__main__":
    main(sys.argv[1:])
    

