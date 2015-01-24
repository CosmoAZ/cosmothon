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
#import time

class lensingClass(object):
    """
    cosmocalcs.cosmologyCalculator(h, 1, .05, .8)
    """

    """
    Initialization of z value by class call
    """
    def __init__(self, h=0.7, omegamat=0.2, omegaDE=0.8):

        self.numbins = 100                # Number of bins to sort redshifts in
	self.integrationbins = 50         # Number of integration steps to determine gaussian probability
	self.z_0 = 0.05               

	self.units = 1                    # 1 for Mpc unit calculations
	self.numgalaxies = 0

	self.h = h                        # Hubble parameters

        self.omegams = [omegamat, omegaDE, 0.05]
        self.wX=-1.                       # dark energy equation of state today
        self.wXa=0.                       # dark energy equation of state evolution
	self.width = .05                  # Redshift width for Gaussian probability function

    """
	Extracts redshifts from DLS data files
    """
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

	# Extract redshift field from file
        #redshiftListF1 = F1datacube[3]
        redshiftListF2 = F2datacube.field(1)
        #redshiftListF3 = F3datacube[3]
        #redshiftListF4 = F4datacube[3]
        #redshiftListF5 = F5datacube[3]

        print "Printing redshift list "
        print redshiftListF2
        print "\n\n"

	self.numgalaxies = len(redshiftListF2)
               
        return redshiftListF2

    """
    Returns the inputted redshift values sorted from smallest to largest.
    """
    def sort_Redshift_Values(self, zlist):

        sortedzlist = sorted(zlist)

        return sortedzlist

    """
    Integrand for weight function integration to be used in lensingWeight Function

    \int_{z}^{\inf} \ frac{D(z, z^')}{D(z^')} n_i(z^') dz^'

    """
    def lensingWeightFunctionIntegrand(self, zprime, z, numList, orgList, zlist):

	"""
	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)
        redShiftCosmocalcs.setEmissionRedShift(zprime)

	"""
	Extraction of angular diameter distance.	
	"""
	if self.units == 1:
	        D_A = redShiftCosmocalcs.AngularDiameterDistanceMpc()
	else:
		D_A = redShiftCosmocalcs.AngularDiameterDistance()

	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)
	redShiftCosmocalcs.setEmissionRedShift(zprime)
	D_A_source = redShiftCosmocalcs.AngularDiameterDistance()

	"""
	Returns integration integrand for lensing weight function calculation
	"""
        #return self.doubletAngularDiameterDistance(z, zprime)/D_A*n_save
	return (D_A*(1+zprime) - D_A_source*(1+z))/(D_A*(1+zprime))*n_save


    """
    Returns integration integrand list.

    \int_{z}^{\inf} \ frac{D(z, z^')}{D(z^')} n_i(z^') dz^'
    """
    def calcIntegrationArray(self, z, orgList, numList, zlist):

	"""
	Calculates array used as integrand in the integration to determine the lensing weight function

        Requires input of the redshift to calculate the weight function with respect to, and the binned 
		redshift values and density list. 
        Outputs the integrand list.
	"""

        outarray = []

	"""
    	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

	# Loops through the binned redshift list, appends to the integration array for redshifts >= the
	#     source redshift
        for zprime in orgList: 

            if zprime >= z:

		"""
		Sets the emission redshift using cosmocalcs class file, will extract angular diameter distance from this.
		"""
	        redShiftCosmocalcs.setEmissionRedShift(zprime)
		D_A = redShiftCosmocalcs.AngularDiameterDistanceMpc()
	
		redShiftCosmocalcs.setEmissionRedShift(z)
		D_A_source = redShiftCosmocalcs.AngularDiameterDistanceMpc()
	
		"""
		Returns integration integrand for lensing weight function calculation
		"""

                outarray.append((D_A*(1+zprime) - D_A_source*(1+z))/(D_A*(1+zprime))*self.n_i(zprime, numList, orgList, zlist)) 
	
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
       	DM2 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistanceMpc()
        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z1)
	DM1 = doubletAngularDiameterDistanceCosmocalcs.TransComovDistanceMpc()	

        """
        Calculates total angular diameter distance between two redshifts
        """
        DA12 = 1/(1 + z2) * (DM2*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5) - DM1*pow(1 + omegamat*pow(DM2, 2)/pow(DH, 2), 0.5))

        return DA12


    """
    Calculates the lensing weight function, W, to be used in the power spectrum calculation
    
    Inputs a set of redshifts, and uses cosmocalcs to calculate the weight at each z
    """
    def lensingWeightFunction(self, z, lensingzlist, orgList, numList):

        cosmocalcsClass = cosmocalcs.cosmologyCalculator(self.h, self.omegams[1], self.omegams[2], self.wX, self.wXa)

        cosmocalcsClass.setEmissionRedShift(z)
        D_A_source = cosmocalcsClass.AngularDiameterDistanceMpc()

	#print "here"
	nlist = []
	for zs in orgList:
		#print z
		nlist.append(self.n_i(zs, numList, orgList, lensingzlist))
	#print nlist
	#sys.exit()
	#print "not here"

        integrand = []
	i = 1
	

	# Loops through the binned redshift list, appends to the integration array for redshifts >= the
	#     source redshift
        for zprime in orgList[1:len(orgList) - 1]: 
            #cosmocalcsClass.setEmissionRedShift(zprime)
            #D_A_source = cosmocalcsClass.AngularDiameterDistanceMpc()
	    #print "Determining lensing weight function at zprime: ", zprime	
	    #print str(zprime) + "\t" + str(z) 
            if zprime >= z:
		#print "MFW"
		"""
		Sets the emission redshift using cosmocalcs class file, will extract angular diameter distance from this.
		"""
	        cosmocalcsClass.setEmissionRedShift(zprime)
		D_A = cosmocalcsClass.AngularDiameterDistanceMpc()
	
		"""
		Returns integration integrand for lensing weight function calculation
		"""
		#print nlist[i]
                integrand.append((D_A*(1+zprime) - D_A_source*(1+z))/(D_A*(1+zprime))*0.68*nlist[i])#*self.n_i(zprime, numList, orgList, lensingzlist)) 

	    i += 1

	x = [i for i in orgList[0:len(orgList)] if i >= z]

	#print integrand

	W = 0

	for i in range(len(integrand) - 1):

		W += (integrand[i+1] + integrand[i])/2.0*(x[i+1] - x[i])

        cosmocalcsClass.setEmissionRedShift(z)
       	D_A = cosmocalcsClass.AngularDiameterDistanceMpc()

        #W *= 3/2*omegamat*1/self.h*(1 + zlim)*D_A*W/3*pow(10, -4)
        W *= 3/2*self.omegams[0]*100*100*self.h*self.h*D_A*W/300000.
        #print "W is ", W#, "\ntime is: ", stuff

        return W 

    """
	Organizes the redshifts into a set of binned arrays.

 	Outputs orgList, a set of binned redshifts, and numList, the number of source galaxies within each bin
    """
    def organizeRedshifts(self, zlist):

	if zlist == []:

		return [], []

	#from numpy import roll

	#print "Number of galaxies found is: ", len(zlist)

        orgList = []
        numList = []
	
        minList = 0.
        #maxList = max(zlist)
	maxList = 5.

        delta = (maxList - minList)/(self.numbins-1)

	binVal = 0

	for i in range(self.numbins):

		orgList.append(binVal)

		#print "Printing binVal ", binVal

		binVal += delta        

		numList.append(0)

	count = 0

	for i in range(len(zlist)):

		for j in range(self.numbins - 1):

			if zlist[i] >= orgList[j] and zlist[i] < orgList[j+1]:

				#if j > 50:
				numList[j] += 1
				count += 1
        return orgList, numList

    """
	Calculates the probability of finding a galaxy of certain redshift within a bin
    """
    def n_i(self, z, numList, orgList, zlist):
	#print "max(orgList): ", max(orgList)
        from numpy import linspace

	nlist = []

	n_z = 0
	deltax = []
	integrand = []
	savedi = 0
	galaxynum = 0

	#print orgList

	# Runs through all the binned redshifts, determines which bin the redshift is associated with
        for i in range(len(orgList) - 1):

            if z >= orgList[i] and z < orgList[i+1]:

		savedi = i
		break
	#print "saved i: " + str(i)
	"""delta = orgList[2] - orgList[1]

	for zt in zlist:
	
		#print zt
	
		if zt >= orgList[savedi] - delta and zt < orgList[savedi+1] + delta:

			#galaxynum += 1

			integrand = []

			x = linspace(orgList[savedi], orgList[savedi], 20)
			
			#x = linspace(min(orgList), max(orgList), 1)

			# Creates gaussian probability array for integration
			for x_prime in x:

				integrand.append(self.gaussFit(x_prime, zt))

			# Integrates this to determine P(z'|z)
			for tempxval in range(len(integrand) - 1):

				n_z += (integrand[tempxval+1] + integrand[tempxval])/2.0*(x[1] - x[0])

			#nlist.append(numList[savedi]/(orgList[1] - orgList[0])*n_z/265601.0)

	f = open('n_z.txt', 'a')
	f.write(str(n_z) + '\n')
	f.close()"""

        #return numList[savedi]/265600.*n_z
	#print "stuff: " + str(numList[savedi]/(orgList[1] - orgList[0])/265601.0/max(orgList))
	return numList[savedi]/(orgList[1] - orgList[0])/265601.0/max(orgList)#self.numgalaxies
	#return 300000*pow(z, 2)*math.exp(-pow(z/.35, 1.2))/14248.9*n_z/self.numgalaxies

    """
	Returns a gaussian probability cetered at z_0 at point x with a given width (normalized)
    """
    def gaussFit(self, x, z_0):

        return math.exp(-pow(x - z_0, 2.)/(2.*pow(self.width, 2)))/(pow(2.*math.pi, 0.5)*self.width)

    def plotWZ(self, Zarray, Warray):

	#Warray =[0.001479565,0.002934595,0.004365565,0.00577293,0.007157145,0.008518655,0.009857885,0.011175265,0.012471195,0.01374609,0.015000335,0.016234315,0.01744841,0.018642985,0.0198184,0.02097501,0.02211315,0.02323317,0.024335385,0.025420125,0.0264877,0.027538425,0.028572595,0.029590505,0.03059245,0.03157871,0.032549555,0.033505265,0.0344461,0.03537232,0.03628418,0.03718193,0.038065815,0.03893607,0.039792925,0.04063662,0.041467375,0.042285405,0.04309093,0.04388416,0.0446653,0.045434555,0.046192125,0.0469382,0.047672975,0.04839664,0.04910937,0.04981135,0.05050275,0.051183755,0.051854525,0.05251523,0.05316603,0.053807085,0.05443856,0.055060595,0.055673355,0.056276975,0.05648664,0.056680845,0.05660439,0.056510525,0.05631738,0.05610734,0.05644715,0.05677662,0.05640809,0.056024255,0.056197255,0.05635945,0.056114835,0.055857225,0.05576731,0.05566597,0.055341955,0.05500652,0.054874285,0.054731915,0.05446744,0.05419311,0.054107225,0.054012655,0.05260803,0.05120618,0.05114498,0.051076545,0.05060548,0.05012755,0.049990985,0.04984783,0.04951141,0.04916875,0.04896569,0.04875673,0.04816497,0.047569095,0.04744164,0.04730912,0.046775435,0.04623811,0.04606945,0.045896295,0.04553726,0.045174345,0.04486793,0.044557755,0.043940745,0.04332237,0.042431565,0.041543745,0.04091384,0.04028372,0.03997518,0.03966411,0.039176245,0.03868732,0.038224945,0.03776157,0.037372175,0.03698144,0.03661773,0.0362527,0.03576438,0.03527608,0.03485722,0.03443793,0.033898025,0.033359305,0.032882105,0.03240558,0.031892455,0.03138073,0.03100628,0.030631835,0.03016733,0.02970405,0.02929534,0.028887385,0.02846382,0.028041385,0.02772194,0.02740263,0.02692371,0.026447105,0.026219035,0.025990595,0.025721835,0.025453105,0.02492072,0.024392215,0.024098585,0.023805475,0.02350273,0.02320071,0.0227684,0.02233881,0.02198747,0.02163776,0.0214193,0.02120103,0.02095359,0.020706695,0.02043721,0.020168595,0.01975409,0.019342855,0.01887996,0.01842164,0.0181032,0.017786725,0.0174159,0.017048155,0.016577715,0.016112965,0.01595517,0.01579764,0.015416475,0.015039255,0.01482707,0.01461588,0.01439162,0.014168595,0.013870875,0.01357578,0.013331355,0.013088685,0.012793195,0.012500605,0.012237615,0.011976995,0.011689945,0.011405955,0.011176695,0.01094941,0.010726235,0.010505025,0.010260235,0.010018,0.009784645,0.00955374,0.00940691,0.009260985,0.009094465,0.008929225,0.00880896,0.00868932,0.008477695,0.008268445,0.00809789,0.007928915,0.007834595,0.007740695,0.007669395,0.0075983,0.00744932,0.00730166,0.007193365,0.00708575,0.006991145,0.00689706,0.00684261,0.006788275,0.006643205,0.006499575,0.00640384,0.006308715,0.00625453,0.006200495,0.006152095,0.006103805,0.006044315,0.00598504,0.005865705,0.00574749,0.00568393,0.00562066,0.00556674,0.005513015,0.00541518,0.005318155,0.005263635,0.005209335,0.005156585,0.00510405,0.00503006,0.00495656,0.00488748,0.004818835,0.00471289,0.00460807,0.004464485,0.004323125,0.004276955,0.004231,0.00416819,0.00410581,0.004021995,0.00393901,0.00386178,0.003785285,0.00371666,0.003648635,0.003593805,0.00353936,0.003473115,0.00340747,0.00335426,0.00330145,0.003260765,0.00322031,0.00314631,0.003073155,0.00302126,0.00296979,0.002936085,0.002902555,0.00283841,0.00277497,0.002714625,0.002654935,0.002581645,0.00250937,0.002471965,0.002434825,0.00240381,0.002372985,0.002213785,0.002060125,0.002031895,0.002003855,0.001979015,0.00195432,0.00194675,0.001939185,0.001908455,0.001877965,0.001866795,0.001855655,0.00178809,0.001721785,0.001564345,0.001414475,0.00139968,0.001384955,0.00137481,0.001364695,0.00129197,0.00122124,0.001201845,0.0011826,0.001179585,0.00117657,0.00108963,0.00100604,0.000986005,0.00096617,0.00092059,0.00087612,0.000872955,0.0008698,0.000866005,0.000862215,0.000860915,0.00085961,0.00080967,0.000761235,0.00075885,0.000756465,0.0007533,0.00075014,0.000728615,0.000707405,0.00068096,0.00065502,0.00065173,0.000648445,0.000627345,0.0006066,0.000590525,0.00057467,0.00056683,0.000559045,0.000537765,0.00051691,0.00051482,0.00051274,0.000494031,0.000475676,0.000473181,0.000470693,0.000448197,0.000426258,0.000376748,0.000330308,0.000310976,0.000292234,0.000281704,0.000271369,0.000263552,0.000255852,0.000222397,4.12955E-05,0.000189171,0.000187058,0.000159928,0.000134931,0.000134194,0.00013346,0.000128975,0.000124567,0.000103312,8.40505E-05,0.000081971,7.99185E-05,0.000079502,7.90865E-05,7.85215E-05,7.79585E-05,7.77665E-05,0.000077575,0.000077274,7.69735E-05,0.000074284,7.16435E-05,0.00007147,7.12965E-05,6.92565E-05,0.000067247,0.000063319,5.95105E-05,5.89715E-05,5.84355E-05,0.000055765,5.31575E-05,4.88148E-05,4.46587E-05,4.44207E-05,4.41833E-05,3.91356E-05,3.4396E-05,3.28316E-05,3.13043E-05,3.12335E-05,3.11628E-05,3.11299E-05,3.10969E-05,3.02286E-05,0.000029373,2.93291E-05,2.92852E-05,2.92291E-05,0.000029173,2.91352E-05,2.90974E-05,2.90778E-05,2.90581E-05,2.90203E-05,2.89825E-05,2.89326E-05,2.88827E-05,2.88388E-05,0.000028795,2.87632E-05,2.87315E-05,2.86697E-05,2.86079E-05,2.85642E-05,2.85205E-05,2.85007E-05,2.84809E-05,2.84253E-05,2.83697E-05,2.83439E-05,2.83182E-05,2.81795E-05,2.80411E-05,2.78854E-05,2.77302E-05,2.77048E-05,2.76793E-05,2.76538E-05,2.76283E-05,2.75969E-05,2.75656E-05,2.74348E-05,2.73043E-05,2.71975E-05,2.70908E-05,2.68168E-05,2.65443E-05,2.55812E-05,2.46363E-05,2.46072E-05,2.45781E-05,2.45215E-05,2.44649E-05,2.42603E-05,2.40566E-05,2.29289E-05,2.18289E-05,2.14713E-05,2.11168E-05,2.10702E-05,2.10236E-05,2.08807E-05,2.07382E-05,2.06669E-05,2.05957E-05,2.05749E-05,2.05541E-05,2.05032E-05,2.04524E-05,2.02119E-05,1.9973E-05,1.98538E-05,0.000019735,1.97197E-05,1.97044E-05,1.96646E-05,1.96248E-05,1.94873E-05,1.93504E-05,1.93353E-05,1.93201E-05,0.000019305,1.92899E-05,1.92457E-05,1.92016E-05,1.9172E-05,1.91424E-05,1.91322E-05,1.91219E-05,1.90828E-05,1.90437E-05,1.90094E-05,1.89752E-05,1.89602E-05,1.89452E-05,1.8935E-05,1.89248E-05,1.89145E-05,1.89043E-05,0.000018894,1.88838E-05,1.88688E-05,1.88538E-05,1.88388E-05,1.88237E-05,1.88135E-05,1.88032E-05,1.8793E-05,1.87827E-05,1.87725E-05,1.87622E-05,1.87519E-05,1.87417E-05,1.87266E-05,1.87116E-05,1.87014E-05,1.86911E-05,1.86808E-05,1.86705E-05,1.86603E-05,1.865E-05,1.86397E-05,1.86294E-05,1.86191E-05,1.86089E-05,1.85986E-05,1.85883E-05,1.8578E-05,1.85677E-05,1.85574E-05,1.85471E-05,1.85368E-05,1.85265E-05,1.85163E-05,1.8506E-05,1.84957E-05,1.84854E-05,1.84751E-05,1.84648E-05,1.84545E-05,1.84442E-05,1.84339E-05,1.84236E-05,1.84134E-05,1.84031E-05,1.83928E-05,1.83825E-05,1.83722E-05,1.83619E-05,1.83517E-05,1.83414E-05,1.83311E-05,1.83208E-05,1.83105E-05,1.83003E-05,1.829E-05,1.82797E-05,1.82694E-05,1.82592E-05,1.82489E-05,1.82386E-05,1.82284E-05,1.82181E-05,1.82079E-05,1.81976E-05,1.81873E-05,1.81771E-05,1.81668E-05,1.81566E-05,1.81464E-05,1.81361E-05,1.81213E-05,1.81065E-05,1.80962E-05,0.000018086,1.80758E-05,1.80656E-05,1.80554E-05,1.80451E-05,1.80304E-05,1.80156E-05,1.80008E-05,1.79861E-05,1.79577E-05,1.79294E-05,1.79101E-05,1.78909E-05,1.78807E-05,1.78706E-05,1.78605E-05,1.78503E-05,1.78266E-05,1.7803E-05,1.77929E-05,1.77828E-05,1.77682E-05,1.77536E-05,1.76445E-05,1.75359E-05,0.000017508,1.74802E-05,1.74703E-05,1.74603E-05,1.73836E-05,0.000017307,1.72883E-05,1.72696E-05,1.72509E-05,1.72323E-05,1.71694E-05,1.71068E-05,1.70882E-05,1.70697E-05,1.69634E-05,1.68575E-05,1.68435E-05,1.68296E-05,1.67851E-05,1.67407E-05,1.66661E-05,1.65916E-05,1.65692E-05,1.65468E-05,1.65115E-05,1.64762E-05,1.64453E-05,1.64145E-05,0.000016345,1.62758E-05,1.62068E-05,1.61379E-05,1.60989E-05,1.60601E-05,1.59535E-05,1.58474E-05,1.583E-05,1.58125E-05,1.57699E-05,1.57274E-05,1.55513E-05,1.53763E-05,1.53468E-05,1.53174E-05,1.51974E-05,1.50779E-05,1.4959E-05,1.48406E-05,1.47875E-05,1.47345E-05,1.44886E-05,1.42448E-05,1.41415E-05,1.40386E-05,1.38185E-05,1.36003E-05,1.35267E-05,1.34534E-05,0.000013338,1.32232E-05,1.2927E-05,1.26343E-05,1.25451E-05,1.24564E-05,1.22061E-05,1.19585E-05,0.000011872,1.17859E-05,1.17253E-05,1.16649E-05,0.000011247,1.0837E-05,1.07448E-05,1.06531E-05,1.04533E-05,1.02555E-05,1.01395E-05,1.00242E-05,9.94255E-06,9.86125E-06,9.6177E-06,9.3773E-06,9.18435E-06,8.9935E-06,8.80475E-06,8.6181E-06,8.50645E-06,8.39555E-06,8.27345E-06,8.1523E-06,8.05575E-06,7.95985E-06,7.7685E-06,7.57955E-06,7.2327E-06,6.89415E-06,6.82215E-06,6.7506E-06,6.65795E-06,6.56605E-06,6.36165E-06,6.16055E-06,6.04425E-06,5.9291E-06,5.8452E-06,5.76195E-06,5.6128E-06,5.46565E-06,5.3661E-06,5.26755E-06,5.13465E-06,5.0035E-06,4.91311E-06,4.82357E-06,4.72815E-06,4.63373E-06,4.50514E-06,4.37842E-06,4.25998E-06,4.14323E-06,4.05516E-06,3.96809E-06,3.84953E-06,3.73283E-06,3.67514E-06,3.61792E-06,3.47415E-06,3.33338E-06,3.25482E-06,3.17725E-06,3.00348E-06,2.83469E-06,2.7951E-06,2.75582E-06,2.6999E-06,2.64458E-06,2.59649E-06,2.54887E-06,2.40832E-06,2.27184E-06,2.1741E-06,2.07857E-06,2.03178E-06,1.98556E-06,1.92136E-06,1.85825E-06,1.84605E-06,1.83389E-06,1.74514E-06,1.65864E-06,1.63537E-06,1.61227E-06,1.60228E-06,1.59233E-06,1.58242E-06,1.57254E-06,1.48445E-06,1.39895E-06,1.38371E-06,1.36855E-06,1.34756E-06,1.32673E-06,1.31543E-06,1.30419E-06,1.29764E-06,1.29112E-06,1.22522E-06,1.16108E-06,1.12664E-06,1.09275E-06,1.09001E-06,1.08728E-06,1.07503E-06,1.06286E-06,1.03105E-06,9.99725E-07,9.49985E-07,9.0154E-07,8.9622E-07,8.9092E-07,8.74205E-07,8.5766E-07,8.4968E-07,8.41745E-07,8.22765E-07,8.04015E-07,7.24705E-07,6.4955E-07,6.4593E-07,6.4232E-07,6.33065E-07,6.23885E-07,6.15575E-07,6.0732E-07,6.0148E-07,5.9567E-07,5.44955E-07,4.96525E-07,4.89863E-07,4.8325E-07,4.8088E-07,4.78516E-07,4.7616E-07,4.7381E-07,4.3613E-07,4.0003E-07,3.94091E-07,3.882E-07,3.74914E-07,3.61867E-07,3.40187E-07,3.19189E-07,3.1731E-07,3.15437E-07,2.98542E-07,2.82121E-07,2.74012E-07,2.66026E-07,2.4705E-07,2.28785E-07,2.26264E-07,2.23758E-07,1.8988E-07,1.58797E-07,1.5277E-07,1.46863E-07,1.39948E-07,1.33203E-07,1.31671E-07,1.30149E-07,1.28996E-07,1.27849E-07,1.15219E-07,1.03252E-07,1.00957E-07,9.86895E-08,9.86355E-08,9.8582E-08,9.7589E-08,9.66015E-08,9.07405E-08,8.5066E-08,8.50195E-08,8.49735E-08,3.3191E-08,8.14285E-08,7.4154E-08,6.72235E-08,6.64135E-08,6.56095E-08,6.50645E-08,6.4522E-08,6.248E-08,6.04715E-08,5.72995E-08,5.4215E-08,5.2805E-08,5.14145E-08,5.13865E-08,5.1359E-08,5.02115E-08,4.90779E-08,4.81755E-08,4.72819E-08,4.70408E-08,4.68004E-08,4.6347E-08,4.58961E-08,4.52364E-08,4.45819E-08,4.16773E-08,3.88721E-08,3.7305E-08,3.57709E-08,3.57517E-08,3.57324E-08,3.55265E-08,3.53212E-08,3.43793E-08,3.34507E-08,3.25352E-08,3.1633E-08,2.95434E-08,2.75263E-08,2.63772E-08,2.52531E-08,2.32428E-08,2.13168E-08,2.07335E-08,2.01585E-08,1.98689E-08,1.95815E-08,1.86181E-08,1.76796E-08,1.74093E-08,1.71413E-08,1.67479E-08,1.63592E-08,1.59754E-08,1.55963E-08,1.38012E-08,1.21168E-08,1.14702E-08,1.08418E-08,1.0836E-08,1.08302E-08,9.82735E-09,8.87375E-09,8.5033E-09,8.141E-09,8.1367E-09,8.13235E-09,7.18525E-09,6.29725E-09,5.83635E-09,5.39315E-09,4.9677E-09,4.55994E-09,4.29715E-09,4.04229E-09,3.9168E-09,3.79335E-09,3.61295E-09,3.43705E-09,2.99247E-09,2.57888E-09,2.43096E-09,2.87485E-10,2.23986E-09,2.19276E-09,2.10124E-09,2.01172E-09,1.96719E-09,1.92318E-09,1.71452E-09,1.51795E-09,1.47946E-09,1.44149E-09,1.44073E-09,1.43997E-09,1.40255E-09,1.36564E-09,1.22503E-09,1.09212E-09,1.05968E-09,1.02775E-09,9.9632E-10,9.6539E-10,8.7655E-10,7.9204E-10,7.3797E-10,6.8584E-10,6.11435E-10,5.41335E-10,5.4105E-10,5.4077E-10,5.182E-10,4.96127E-10,4.1338E-10,3.38214E-10,3.38038E-10,3.37862E-10,3.03076E-10,2.70196E-10,2.39219E-10,2.10142E-10,1.82962E-10,1.57676E-10,1.34282E-10,1.12776E-10,1.02703E-10,9.3106E-11,9.3058E-11,9.30095E-11,6.7165E-11,4.55277E-11,4.55041E-11,4.54806E-11,3.33971E-11,2.31805E-11,1.87665E-11,1.48202E-11,1.48126E-11,1.48049E-11,1.13292E-11,8.3192E-12,5.77425E-12,3.69361E-12,2.07659E-12,9.2245E-13,2.30495E-13,0,0,0,0,0]

	#from numpy import linspace
	#Zarray = linspace(0, max(Zarray), len(Warray))

	import matplotlib.pyplot as plt

	plt.plot(Zarray, Warray, 'g--')
	plt.xlabel('Redshift')
   	plt.ylabel('W(z)')
    	plt.title('Lensing Weight Function vs. Redshift')
    	plt.axis([min(Zarray), max(Zarray), min(Warray), max(Warray)])
    	plt.grid(True)
    	plt.savefig('WvsZ.png')
    	plt.close()

    def n_list(self, numList, orgList, zlist):

        from numpy import linspace
	delta = orgList[2] - orgList[1]
	nlist = []
	n_z = 0

	for z in orgList:

            for i in range(len(orgList) - 1):

	        if z >= orgList[i] and z < orgList[i+1]:

		    savedi = i
		    break

	    for zt in zlist:

		if zt >= orgList[savedi] - delta and zt < orgList[savedi+1] + delta:

			integrand = []

			x = linspace(orgList[savedi], orgList[savedi+1], 100)
			#x = linspace(min(orgList), max(orgList), 10000)

			# Creates gaussian probability array for integration
			for x_prime in x:

				integrand.append(self.gaussFit(x_prime, zt))

			# Integrates this to determine P(z'|z)
			for tempxval in range(len(integrand) - 1):

				n_z += (integrand[tempxval+1] + integrand[tempxval])/2.0*(x[tempxval+1] - x[tempxval])

	    nlist.append(numList[savedi]/(orgList[1] - orgList[0])/0.68*265601.0)#*n_z)

	return nlist

"""



        from numpy import linspace

	nlist = []

	n_z = 0
	deltax = []
	integrand = []
	savedi = 0
	galaxynum = 0

	# Runs through all the binned redshifts, determines which bin the redshift is associated with
        for i in range(len(orgList) - 1):

            if z >= orgList[i] and z < orgList[i+1]:

		savedi = i
		break

	delta = orgList[2] - orgList[1]

	for zt in zlist:
	
		if zt >= orgList[savedi] - delta and zt < orgList[savedi+1] + delta:

			galaxynum += 1

			integrand = []

			x = linspace(orgList[savedi], orgList[savedi+1], 100)
			#x = linspace(min(orgList), max(orgList), 10000)

			# Creates gaussian probability array for integration
			for x_prime in x:

				integrand.append(self.gaussFit(x_prime, zt))

			# Integrates this to determine P(z'|z)
			for tempxval in range(len(integrand) - 1):

				n_z += (integrand[tempxval+1] + integrand[tempxval])/2.0*(x[tempxval+1] - x[tempxval])

			nlist.append(numList[savedi]/(orgList[1] - orgList[0])*n_z/265601.0)

        #return numList[savedi]/265600.*n_z
	#return numList[savedi]/(orgList[1] - orgList[0])*n_z/265601.0#/max(orgList)#self.numgalaxies
	#return 300000*pow(z, 2)*math.exp(-pow(z/.35, 1.2))/14248.9*n_z/self.numgalaxies
	return nlist"""

