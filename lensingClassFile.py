## @package lensingClass
# Contains lensing calculation functions
#
# lensing class

import sys
import cosmocalcs
import pyfits
import copy
#import numdisplay

class lensingClass(object):
    """
    cosmocalcs.cosmologyCalculator(h, 1, .05, .8)
    """

    """
    Initialization of z value by class call
    """
    def __init__(self, z):

        self.z = z
        self.numbins = 100

    def getRedshiftsFromDLS(self):
        
        # wcs info
        """fileF1 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F1wcs.fits')
        hdrF1 = fileF1[0].header
        fileF2 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F2wcs.fits')
        hdrF2 = fileF2[0].header
        fileF3 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F3wcs.fits')
        hdrF3 = fileF3[0].header
        fileF4 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F4wcs.fits')
        hdrF4 = fileF4[0].header
        fileF5 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F5wcs.fits')
        hdrF5 = fileF5[0].header"""
        
        #F1 = pyfits.open('/home/akilgall/field1_sources.fits')
        F2 = pyfits.open('/home/akilgall/field2_sources.fits')
        #F3 = pyfits.open('/home/akilgall/field3_sources.fits')
        #F4 = pyfits.open('/home/akilgall/field4_sources.fits')
        #F5 = pyfits.open('/home/akilgall/field5_sources.fits')
        
        #F1header = pyfits.getheader('/home/akilgall/field1_sources.fits')
        F2header = pyfits.getheader('/home/akilgall/field2_sources.fits')
        #F3header = pyfits.getheader('/home/akilgall/field3_sources.fits')
        #F4header = pyfits.getheader('/home/akilgall/field4_sources.fits')
        #F5header = pyfits.getheader('/home/akilgall/field5_sources.fits')
        
        #F1datacube = pyfits.getdata('/home/akilgall/field1_sources.fits')
        F2datacube, F2header = pyfits.getdata('/home/akilgall/field2_sources.fits', header=True)
        #F3datacube = pyfits.getdata('/home/akilgall/field3_sources.fits')
        #F4datacube = pyfits.getdata('/home/akilgall/field4_sources.fits')
        #F5datacube = pyfits.getdata('/home/akilgall/field5_sources.fits')

        print "Printing table headers"
        print F2header
        print "\n\n"

        print "Printing data cubes"
        print F2datacube
        print "\n\n"

        #redshiftListF1 = F1datacube[3]
        redshiftListF2 = F2datacube.field(2)
        #redshiftListF3 = F3datacube[3]
        #redshiftListF4 = F4datacube[3]
        #redshiftListF5 = F5datacube[3]

        print "Printing redshift list "
        print redshiftListF2
        print "\n\n"
        
       #numdisplay.display(redShiftListF2)
                
        return redshiftListF2



    def sort_Redshift_Values(self, zlist):

        sortedzlist = sorted(zlist)

        return sortedzlist

    """
    Integrand for spectrum function integration to be used in scipy quad function
    """
    def lensingWeightFunctionIntegrand(self, x, z, numList, orgList):

        omegams = [1, 0.05, 0.2]
        wX=-1. # dark energy equation of state today
        wXa=0. # dark energy equation of state evolution
        h = 0.7

	"""
	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(h, omegams[1], omegams[2], wX, wXa)

	"""
	Sets the emission redshift using cosmocalcs class file, will extract angular diameter distance from this.
	"""
        redShiftCosmocalcs.setEmissionRedShift(z)

	"""
	Extraction of angular diameter distance.	
	"""
        D_A = redShiftCosmocalcs.AngularDiameterDistance()

	"""
	Returns integration integrand for power spectrum calculation
	"""
        return self.doubletAngularDiameterDistance(z, x)/D_A*self.n_i(z, numList, orgList)

    def calcIntegrationArray(self, init_cosmologyCalculator, z, zlist, numList, orgList):

        outarray = []

        for zval in zlist:

            if zval > z:

                print "In calcIntegrationArray, checking z: ", zval

                outarray.append(self.lensingWeightFunctionIntegrand(zval, z, numList, orgList)) 

        return outarray

    """
    Calculates Angular Diameter distance between z1 and z2
    """
    def doubletAngularDiameterDistance(self, z1, z2):

	"""
	Initializes cosmocalcs class
	"""
	doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(1, 1, 0, 0)

        omegamat = 0.7

        """
        Calculates and then returns angular diameter distance for z2
        """
        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z2)
        DM2 = doubletAngularDiameterDistanceCosmocalcs.AngularDiameterDistance()

        doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(1, 1, 0, 0)

        """
        Calculates and then returns angular diameter distnace for z1
        """
        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z1)
        DM1 = doubletAngularDiameterDistanceCosmocalcs.AngularDiameterDistance()

        """
        """
        DH = 1/(1 + z2)

        """
        Calculates total angular diameter distance between two redshifts
        """
        DA12 = 1/(1 + z2) * (DM2*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5) - DM1*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5))

        return DA12

    """
    Calculates the lensing weight function, W, to be used in the power spectrum calculation
    
    Inputs a set of redshifts, and uses cosmocalcs to calculate the weight at each z
    """
    def lensingWeightFunction(self, lensingzlist):

        omegams = [1, 0.05, 0.2] 
        wX=-1. # dark energy equation of state today
        wXa=0. # dark energy equation of state evolution

        print "Beginning of lensingWeightFunction\n"

        from scipy.integrate import simps

        f = open('Wvalues.txt', 'w')

        n_i = 1
        h = 0.7
        omegamat = 0.7
         
        init_cosmologyCalculator = cosmocalcs.cosmologyCalculator(h, omegams[1], omegams[2], wX, wXa)
               
        for z in lensingzlist:

            print "z is: ", z 

            init_cosmologyCalculator.setEmissionRedShift(z)
            D_A = init_cosmologyCalculator.AngularDiameterDistance()

            print "calculating integrand array"

            #integrand = self.calcIntegrationArray(z, lensingzlist)
            integrand = self.calcIntegrationArray(init_cosmologyCalculator, z, lensingzlist)

            print "integrating now"

            W = simps(integrand, [i for i in lensingzlist if i >= z]) 

            W *= 3/2*omegamat*h*h*(1 + self.z)*D_A*W[1]

            print "W is ", W

            f.write(W)
            f.write('\n')

        f.close()

        return W

    def organizeRedshifts(self, zlist):

        orgList = []
        numList = []

        minList = min(zlist)
        maxList = max(zlist)

        delta = (maxList - minList)/self.numbins

        listVal = minList

        for i in range(1, 100):

            orgList.append(listVal)
            numList.append(len([val for val in listVal if val >= listVal and val < listVal + delta]))

            listVal += delta

        return orgList, numList

        #        integrand = lambda x: cosmocalcs.calcAngularDistance(z, x)/cosmocalcs.calcAngularDistance(z)#*n_i(x)

    def n_i(self, z, numList, orgList):

        delta = orgList[1] - orgList[0]

        low = min([i for i in orgList if i >= orgList and i < orgList + delta])

        return low


