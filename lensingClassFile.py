# @package lensingClass
# Contains lensing calculation functions
#
# lensing class

import sys
import cosmocalcs
import pyfits
import copy
import powerspec
import math
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
        self.numbins = 300

	self.h = 0.7

        self.omegams = [0.3, 0.7, 0.05]
        self.wX=-1. # dark energy equation of state today
        self.wXa=0. # dark energy equation of state evolution
#        self.width = .01
        self.width = .01

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
        #F2header = pyfits.getheader('/home/akilgall/field2_sources.fits')
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
        redshiftListF2 = F2datacube.field(3)
        #redshiftListF3 = F3datacube[3]
        #redshiftListF4 = F4datacube[3]
        #redshiftListF5 = F5datacube[3]

        print "Printing redshift list "
        print redshiftListF2
        print "\n\n"
               
        return redshiftListF2

    """
    Returns the inputted redshift values sorted from smallest to largest.
    """
    def sort_Redshift_Values(self, zlist):

        sortedzlist = sorted(zlist)

        return sortedzlist

    """
    Integrand for spectrum function integration to be used in scipy quad function
    """
    def lensingWeightFunctionIntegrand(self, x, z, numList, orgList):

	"""
	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

	"""
	Sets the emission redshift using cosmocalcs class file, will extract angular diameter distance from this.
	"""
        redShiftCosmocalcs.setEmissionRedShift(x)

	"""
	Extraction of angular diameter distance.	
	"""
        D_A = redShiftCosmocalcs.AngularDiameterDistance()

	"""
	Returns integration integrand for power spectrum calculation
	"""
        return self.doubletAngularDiameterDistance(z, x)/D_A*self.n_i(x, numList, orgList)



    """
    Returns integration integrand list.
    """
    def calcIntegrationArray(self, z, orgList, numList):

	"""
	Calculates array used as integrand in the integration to determine the lensing weight function

        Requires input of the redshift to calculate the weight function with respect to, and the binned redshift values and density list. 
        Outputs the integrand list.
	"""

        outarray = []

#	print "Printing orgList: ", orgList

        for zval in orgList: #zlist:

            if zval > z:

 #               print "In calcIntegrationArray, checking z: ", zval

                outarray.append(self.lensingWeightFunctionIntegrand(zval, z, numList, orgList)) 

        return outarray



    """
    Calculates Angular Diameter distance between z1 and z2
    """
    def doubletAngularDiameterDistance(self, z1, z2):

	"""
	Initializes cosmocalcs class
	"""
	doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

        omegamat = 0.3

        """
        Calculates and then returns angular diameter distance for z2
        """
        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z2)
        DM2 = doubletAngularDiameterDistanceCosmocalcs.AngularDiameterDistance()

        doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

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
    def lensingWeightFunction(self, zlim, lensingzlist, orgList, numList):

        print "Beginning of lensingWeightFunction\n"

        from scipy.integrate import simps

        omegamat = 0.3
                        
	#orgList, numList = self.organizeRedshifts(lensingzlist)

        """for z in lensingzlist[1:len(lensingzlist)-1]:

            print "z is: ", z 

            init_cosmologyCalculator = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

            init_cosmologyCalculator.setEmissionRedShift(z)
            D_A = init_cosmologyCalculator.AngularDiameterDistance()

            print "calculating integrand array"

            #integrand = self.calcIntegrationArray(z, lensingzlist)
            integrand = self.calcIntegrationArray(z, lensingzlist, numList)

            print "integrating now"

            print "Size of integrand is ", len(integrand)

            print "Size of x vals is ", len([i for i in orgList[0:len(orgList)-1] if i >= z])

            print "\n\n\n\n"

            print "Integrand is: ", integrand
            print "x is: ", [i for i in orgList if i >= z]

            W = simps(integrand, [i for i in orgList[0:len(orgList)-1] if i >= z]) 

            W *= 3/2*omegamat*self.h*self.h*(1 + self.z)*D_A*W

            print "W is ", W

            f.write(str(z) + "\t" + str(W) + "\n")

        f.close()"""


        print "z is: ", zlim

        init_cosmologyCalculator = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

        init_cosmologyCalculator.setEmissionRedShift(zlim)
        D_A = init_cosmologyCalculator.AngularDiameterDistance()

        print "calculating integrand array"

        integrand = self.calcIntegrationArray(zlim, lensingzlist, numList)

        print "integrating now"

        print "Size of integrand is ", len(integrand)

        print "Size of x vals is ", len([i for i in orgList[0:len(orgList)-1] if i >= zlim])

        print "\n\n\n\n"

#        print "Integrand is: ", integrand
#        print "x is: ", [i for i in orgList if i >= zlim]

        W = simps(integrand, [i for i in orgList[0:len(orgList)-1] if i >= zlim]) 

        W *= 3/2*omegamat*self.h*100*100*self.h*(1 + zlim)*D_A*W/300000

        print "W is ", W

        return W

    def organizeRedshifts(self, zlist):

        orgList = []
        numList = []
	
        minList = min(zlist)
        maxList = max(zlist)

        delta = (maxList - minList)/self.numbins

	binVal = minList

	for i in range(self.numbins):

		orgList.append(binVal)

		print "Printing binVal ", binVal

		binVal += delta        

		numList.append(0)

	for i in range(len(zlist)):

		for j in range(self.numbins - 1):

			if zlist[i] >= orgList[j] and zlist[i] < orgList[j+1]:

				numList[j] += 1

        """for i in range(1, 100):

            orgList.append(listVal)
            #numList.append(len([val for val in listVal if val >= listVal and val < listVal + delta]))

	    for i in range(listVal):

		if 

            #listVal += delta
	"""

#	print "Printing orgList !", orgList

        return orgList, numList

        #        integrand = lambda x: cosmocalcs.calcAngularDistance(z, x)/cosmocalcs.calcAngularDistance(z)#*n_i(x)

    def n_i(self, z, numList, orgList):

        from numpy import linspace
        from scipy.integrate import simps

#        delta = orgList[1] - orgList[0]

#        low = min([i for i in orgList if i >= orgList and i < orgList + delta])

        n_z = -100000000

        for i in range(self.numbins - 1):

            if z >= orgList[i] and z < orgList[i+1]:

                x = linspace(orgList[i], orgList[i+1], 20)

                integrand = []

                for x_prime in x:
                
                    integrand.append(self.gaussFit(x_prime, z))

                n_z = simps(integrand, x)

        return n_z*pow(z, 2)*math.exp(-pow(z/0.25, 2))
#        return pow(z, 2)*math.exp(-pow(z/0.25, 2))


    def gaussFit(self, x, z_0):

        return math.exp(-pow(x - z_0, 2)/(2*pow(self.width, 2)))

    def calcPowerSpectrum(self):


	l = 1

        ps = powerspec.transferFunction(1, 1, 1)

	f = open('PowerSpectrum_results.txt', 'w')

	f.write("z" + '\t' + "Wval" + "\t" + "PS" + '\n')

	lc = lensingClassFile.lensingClass(0)

	F2redshifts = lc.getRedshiftsFromDLS()

	organizedList, numberedList = lc.organizeRedshifts(F2redshifts)

	sortedF2Redshifts = lc.sort_Redshift_Values(organizedList)

	Wval = []

	i = 1

	for z in sortedF2Redshifts:

		print "Determining lensing weight function at z: ", z
        
		Wval.append(lc.lensingWeightFunction(sortedF2Redshifts, organizedList, numberedList))
        
		f.write(z + '\t' + Wval[i] + '\n')

		i += 1      

	f.close()

        return Wval






