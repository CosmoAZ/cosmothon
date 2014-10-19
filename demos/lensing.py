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

    orgList, numList = lc.organizeRedshifts(F2redshifts)

    print "Sorting redshift values"

    sortedF2Redshifts = lc.sort_Redshift_Values(orgList)

    for z in sortedF2Redshifts:

        print "Determining lensing weight function at z: ", z
        
        Wval = lc.lensingWeightFunction(sortedF2Redshifts)
        
        f.write(z + '\t' + Wval + '\n')
      
    f.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    

