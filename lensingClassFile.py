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
        self.numbins = 100
	self.integrationbins = 50
	self.z_0 = 0.05

	self.units = 0 #1 for Mpc

	self.h = 0.7

        self.omegams = [0.3, 0.7, 0.05]
        self.wX=-1. # dark energy equation of state today
        self.wXa=0. # dark energy equation of state evolution
#        self.width = .01
        self.width = .05

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
        redshiftListF2 = F2datacube.field(1)
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
    Integrand for weight function integration to be used in lensingWeight Function
    """
    def lensingWeightFunctionIntegrand(self, zprime, z, numList, orgList):

	"""
	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

	"""
	Sets the emission redshift using cosmocalcs class file, will extract angular diameter distance from this.
	"""
        redShiftCosmocalcs.setEmissionRedShift(zprime)

	"""
	Extraction of angular diameter distance.	
	"""
	if self.units == 1:
	        D_A = redShiftCosmocalcs.AngularDiameterDistanceMpc()
	else:
		D_A = redShiftCosmocalcs.AngularDiameterDistance()

	oi = open('integrand_components.txt', 'a')

        n_save = self.n_i(zprime, numList, orgList)
	
	oi.write('z: ' + str(z) + '\t' + 'z\': ' + str(zprime) + '\t' + 'D_Z: ' + str(D_A) + '\t' + 'DD_A: ' + str(self.doubletAngularDiameterDistance(z, zprime)) +'\t' + 'n_i: ' + str(n_save) + '\n')

	oi.close()

	"""
	Returns integration integrand for power spectrum calculation
	"""
        return self.doubletAngularDiameterDistance(z, zprime)/D_A*n_save



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

        for zprime in orgList: #zlist:

            if zprime >= z:

 #               print "In calcIntegrationArray, checking z: ", zval

                outarray.append(self.lensingWeightFunctionIntegrand(zprime, z, numList, orgList)) 

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

	if self.units == 1:
        	DM2 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistanceMpc()
	else:
		DM2 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistance()

        doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

        """
        Calculates and then returns angular diameter distnace for z1
        """
        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z1)

	if self.units == 1:
        	DM1 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistanceMpc()
	else:
        	DM1 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistance()		

        """
        """

	if self.units == 1:	
        	DH = 100.0/self.h
	else:
		DH = 100.0/self.h

        """
        Calculates total angular diameter distance between two redshifts
        """
        DA12 = 1/(1 + z2) * (DM2*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5) - DM1*pow(1 + omegamat*pow(DM2, 2)/pow(DH, 2), 0.5))

#	f4 = open('dd.txt', 'a')
#	f4.write(str(DA12) + '\n')
#	f4.close()

        return DA12


    """
    Calculates the lensing weight function, W, to be used in the power spectrum calculation
    
    Inputs a set of redshifts, and uses cosmocalcs to calculate the weight at each z
    """
    def lensingWeightFunction(self, zlim, lensingzlist, orgList, numList):

        print "Beginning of lensingWeightFunction\n"

        #from scipy.integrate import simps

        omegamat = 0.3

        print "z is: ", zlim

        init_cosmologyCalculator = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

        init_cosmologyCalculator.setEmissionRedShift(zlim)

	if self.units == 1:
        	D_A = init_cosmologyCalculator.AngularDiameterDistanceMpc()
	else:
        	D_A = init_cosmologyCalculator.AngularDiameterDistance()

        print "calculating integrand array"

        integrand = self.calcIntegrationArray(zlim, lensingzlist, numList)

        print "integrating now"
        print "Size of integrand is ", len(integrand)
        print "Size of x vals is ", len([i for i in orgList[0:len(orgList)-1] if i >= zlim])
        print "\n\n\n\n"

	x = [i for i in orgList[0:len(orgList)] if i >= zlim]

	xfile = open('x.txt', 'a')
	for xxt in x:
		xfile.write(str(xxt) + '\t')
	xfile.write('\n')
	xfile.close()

	W = 0

	for i in range(len(integrand) - 1):

		W += (integrand[i+1] + integrand[i])/2.0*(x[i+1] - x[i])

        #W = simps(integrand, [i for i in orgList[0:len(orgList)-1] if i >= zlim]) 

	f1 = open('W.txt', 'a')
	f1.write('W_i: ' + str(W) + '\t' + 'D_A: ' + str(D_A) + '\t' + '3/2*omegamat*(1+z): ' + str(3/2*omegamat*(1+zlim)) + '\n')
	f1.close()

        #W *= 3/2*omegamat*1/self.h*(1 + zlim)*D_A*W*3*pow(10, -4)
        W *= 3/2*omegamat*(1 + zlim)*D_A*W
        print "W is ", W

        return W

    def organizeRedshifts(self, zlist):

	from numpy import roll

	print "Number of galaxies found is: ", len(zlist)

        orgList = []
        numList = []
	
        minList = min(zlist)
        maxList = max(zlist)

        delta = (maxList - minList)/(self.numbins-1)

	binVal = minList

	for i in range(self.numbins):

		orgList.append(binVal)

		print "Printing binVal ", binVal

		binVal += delta        

		numList.append(0)

	count = 0

	for i in range(len(zlist)):

		for j in range(self.numbins - 1):

			if zlist[i] >= orgList[j] and zlist[i] < orgList[j+1]:

				numList[j] += 1
				count += 1


	ff = open('hist.txt', 'w')

	ff.write('orgList: ' + str(orgList) + '\n')
	ff.write('numList: ' + str(numList) + '\n')
	ff.write('Number of galaxies found: ' + str(len(zlist)) + '\n')
	ff.write('Number of galaxies sorted: ' + str(count) + '\n')

	ff.close()

	#numList = roll(numList, 50)

        return orgList, numList


    def n_i(self, z, numList, orgList):

	f2 = open('n_i.txt', 'a')

        from numpy import linspace

	n_z = 0
	deltax = []
	integrand = []
	savedi = 0

        for i in range(len(orgList) - 1):

            if z >= orgList[i] and z < orgList[i+1]:

		savedi = i

		f2.write('i: ' + str(i) + '\t' + 'orgList[i]: ' + str(orgList[i]) + '\t' + 'orgList[i+1]: ' + str(orgList[i+1]) + '\t' + 'numList[savedi]/265600/(orgList[1] - orgList[0]): ' + str(numList[savedi]/265600.*1/(orgList[1] - orgList[0])) + '\t')

                ##x = linspace(orgList[i], orgList[i+1], 50)
		x = linspace(min(orgList), max(orgList), 10000)

                for x_prime in x:
                
		    if x_prime >= orgList[i] and x_prime <= orgList[i+1]:

		        integrand.append(self.gaussFit(x_prime, z))

		for tempxval in range(len(integrand) - 1):

			n_z += (integrand[tempxval+1] + integrand[tempxval])/2.0*(x[tempxval+1] - x[tempxval])

                #n_z = simps(integrand, x)
		deltax = x[tempxval+1] - x[tempxval];

		break


	f2.write('n: ' + str(numList[savedi]) + '\t' + 'z: ' + str(z) + '\t' + 'n_z: ' + str(n_z) + '\t' + 'deltax: ' + str(deltax) + '\t')

	f2.write('integrand[i]: ')

#	for intg in integrand:

#		f2.write(str(intg) + '\t')

	f2.write('\n\n')
	f2.close()

        #return 4000*n_z*pow(z, 2)*math.exp(-pow(z/self.z_0, 1.2))
	#return numList[savedi]/265600.*1./(orgList[1] - orgList[0])*n_z
	return numList[savedi]/265600.

    def gaussFit(self, x, z_0):

        return math.exp(-pow(x - z_0, 2)/(2*pow(self.width/2, 2)))/(pow(2*math.pi, 0.5)*self.width/2)

    def calcPowerSpectrum(self):


	l = 1

        ps = powerspec.transferFunction(1, 1, 1)

 	f = open('results.txt', 'w')

	print "Initializing lc class\n"

    	lc = lensingClassFile.lensingClass(0)

	print "Extracting redshifts from DLS data"

	F2redshifts = lc.getRedshiftsFromDLS()

	organizedList, numberedList = lc.organizeRedshifts(F2redshifts)

	print "Printing the organized list after organizeRedshifts function call"
	print "OrgList is: ?", organizedList

	print "Sorting redshift values"

	sortedF2Redshifts = lc.sort_Redshift_Values(organizedList)

	for z in sortedF2Redshifts[1:len(sortedF2Redshifts) - 1]:

		print "Determining lensing weight function at z: ", z
        
		Wval = lc.lensingWeightFunction(z, sortedF2Redshifts, organizedList, numberedList)
        
		f.write(str(z) + '\t' + str(Wval) + '\n')
      
	f.close()






