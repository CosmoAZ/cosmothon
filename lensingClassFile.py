## @package lensingClass
# Contains lensing calculation functions
#
# lensing class

import cosmocalcs

class lensingClass(object):

    """
    cosmocalcs.cosmologyCalculator(h, 1, .05, .8)
    """

    print "Before __init__\n"
   
    """
    Initialization of z value by class call
    """
    def __init__(self, z):

        self.z = z

    """
    Integrand for spectrum function integration to be used in scipy quad function
    """
    def integrand(self, x, z):

	i = 1
        print "i"
	i++

	"""
	Initialization of cosmocalcs class for each integration step.  I'm sure that there's a better way to do this.
	"""
	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(1, 1, 0, 0)

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
        return self.doubletAngularDiameterDistance(z, x)/D_A*self.n_I(x)


    """
    Calculates Angular Diameter distance between z1 and z2
    """
    def doubletAngularDiameterDistance(self, z1, z2):

	"""
	Initializes cosmocalcs class
	"""
	doubletAngularDiameterDistanceCosmocalcs = cosmocalcs.cosmologyCalculator(1, 1, 0, 0)

        omegamat = 0.7

        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z2)
        DM2 = doubletAngularDiameterDistanceCosmocalcs.AngularDiameterDistance()

        doubletAngularDiameterDistanceCosmocalcs.setEmissionRedShift(z1)
        DM1 = doubletAngularDiameterDistanceCosmocalcs.AngularDiameterDistance()

        DH = 1/(1 + z2)

        DA12 = 1/(1 + z2) * (DM2*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5) - DM1*pow(1 + omegamat*pow(DM1, 2)/pow(DH, 2), 0.5))

        return DA12


    def lensingWeightFunction(self, z):

        print "Beginning of lensingWeightFunction\n"

        from scipy.integrate import quad

        n_i = 1
        h = 1
        temp_h = h
        omegamat = 0.7

        print "h is " 
        print h

        init_cosmologyCalculator = cosmocalcs.cosmologyCalculator(1, 0.05, 0.08)
        
        init_cosmologyCalculator.setEmissionRedShift(z)
        D_A = init_cosmologyCalculator.AngularDiameterDistance()

        print D_A

        zmax = 50
        zmin = self.z
        b = 0

        W = quad(self.integrand, zmin, zmax, args=(z)) 
	
	W = 3/2*omegamat*h*h*(1 + self.z)*D_A*W[1]

        print "W is ", W

        return W


        #        integrand = lambda x: cosmocalcs.calcAngularDistance(z, x)/cosmocalcs.calcAngularDistance(z)#*n_i(x)


