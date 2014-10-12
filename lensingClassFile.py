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
   
    def __init__(self, z):

        self.z = z

    def integrand(self, x, z):

        print "i"

	redShiftCosmocalcs = cosmocalcs.cosmologyCalculator(1, 1, 0, 0)

        redShiftCosmocalcs.setEmissionRedShift(z)

        D_A = redShiftCosmocalcs.AngularDiameterDistance()

        redShiftCosmocalcs.setEmissionRedShift(z)

        return self.doubletAngularDiameterDistance(z, x)/redShiftCosmocalcs.AngularDiameterDistance()

    def doubletAngularDiameterDistance(self, z1, z2):

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

