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

    print "Determining lensing weight function\n"

    wat = lc.lensingWeightFunction(z)

    print wat

if __name__ == "__main__":
    main(sys.argv[1:])
    

