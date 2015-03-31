import sys
import lensingClassFile
import cosmocalcs
from numpy import linspace
from numpy import logspace
from scipy.interpolate import interp1d
from scipy.integrate import quad
import powerspec
import math

class shearCalcClass(object):

	def __init__(self, shearBins=1, i=0, j=1, h=0.7, omegamat=0.3, omegaDE=0.7, nlistp=[0]):
		self.i = i
		self.j = j
		self.shearBins = shearBins
		self.h = h
		self.omegamat = omegamat
		self.omegaDE = omegaDE
		self.ns = 0.96
		self.DelRsq = 2.1e-9
		self.nummulipoles = 20000
		self.nlistp = nlistp
		#print("\n\nPrinting nlistp\n\n")
		#print nlistp
		#print("\n\nPrinting nlistp\n\n")
		
	def shearPSPrePrep(self, F2redshifts):
		print("\nshearClass:  Omegamat:\t" + str(self.omegamat) + "  OmegaDE:\t" + str(self.omegaDE))
		if self.i > self.shearBins - 1 or self.j > self.shearBins - 1 or self.i < 0 or self.j < 0:
			print "i or j (bin number) is greater than the number of Shear bins: ", self.shearBins
			return

		lc = lensingClassFile.lensingClass(self.h, self.omegamat, self.omegaDE, self.nlistp)

		larray = logspace(math.log10(1), math.log10(1e8), self.nummulipoles)
		pKij = []

		#F2redshifts = lc.getRedshiftsFromDLS()    	
		#organizedList, numberedList = lc.organizeRedshifts(F2redshifts)
		
		#f_x = interp1d(linspace(0, max(zttiarray), len(Warraytt)), Warraytt)

		#redshifts = linspace(0, max(zttiarray), len(organizedList))
		#Warraytt = f_x(redshifts)
		
		Warrayi = []
		Warrayj = []

		shearbinarray = linspace(0, max(F2redshifts), self.shearBins+1)

		if self.shearBins == 1:

			zttiarray = F2redshifts
			organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
			sortedzttiarray = lc.sort_Redshift_Values(organizedListi)

		    	Warrayi = []
			#print "sortedzttiarray: " + sortedzttiarray
		
		    	for z in sortedzttiarray[0:len(sortedzttiarray) - 1]:
				#print "z                                 : " + str(z)
				Wvali = lc.lensingWeightFunction(z, zttiarray, organizedListi, numberedListi)
				Warrayi.append(Wvali)
		else:

			zttiarray = [ii for ii in F2redshifts if ii >= shearbinarray[self.i] and ii < shearbinarray[self.i+1]]
			zttjarray = [ii for ii in F2redshifts if ii >= shearbinarray[self.j] and ii < shearbinarray[self.j+1]]

			organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
			sortedzttiarray = lc.sort_Redshift_Values(organizedListi)			
			
			organizedListj, numberedListj = lc.organizeRedshifts(zttjarray)
			sortedzttjarray = lc.sort_Redshift_Values(organizedListj)

		    	Warrayi = []
		    	Warrayj = []
		
		    	for z in sortedzttiarray[1:len(sortedzttiarray) - 1]:
			
				Wvali = lc.lensingWeightFunction(z, sortedzttiarray, organizedListi, numberedListi)
				Wvalj = lc.lensingWeightFunction(z, sortedzttjarray, organizedListj, numberedListj)

				Warrayi.append(Wvali)
				Warrayj.append(Wvalj)

		#print "Size of Warray: ", len(Warrayi), "\nSize of sortedzttiarray: ", len(sortedzttiarray[0:len(Warrayi)])
		lc.plotWZ(sortedzttiarray[0:len(Warrayi)], Warrayi)

		if (self.omegamat < lc.omegams[2]):	
			self.omegamat = lc.omegams[2]	

		cosmoCalcs = cosmocalcs.cosmologyCalculator(lc.h, self.omegamat, self.omegaDE, lc.wX, lc.wXa)
		tf = powerspec.transferFunction(self.omegamat-lc.omegams[2], lc.omegams[2], self.h)
		gro = powerspec.growthFunctionLCDM(self.omegamat, self.omegaDE)
		ps = powerspec.powerSpectrum(tf, gro, self.ns, self.DelRsq)

		ADDM = []
		for z in sortedzttiarray[1:len(zttiarray) - 1]:
			cosmoCalcs.setEmissionRedShift(z)					
			ADDM.append(cosmoCalcs.AngularDiameterDistanceMpc())

		for l in larray:

			pKij_l = 0
			ii = 0

			for z in sortedzttiarray[1:len(zttiarray) - 1]:
				#print "Determining cross-spectra weight function at z: ", z
				#cosmoCalcs.setEmissionRedShift(z)

				D_AMpc = (1 + z)*ADDM[ii]
				Pk = ps.getPk(l/D_AMpc, z) # in (Mpc/h)^3 units

				if self.shearBins == 1:
					pKij_l += Warrayi[ii]*Warrayi[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegamat*pow(1 + z, 3) + (1 - self.omegaDE - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])
				else:
					pKij_l += Warrayi[ii]*Warrayj[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegaDE*pow(1 + z, 3) + (1 - self.omegaDE - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])

				ii += 1

			pKij.append(pKij_l)

			#f = open('Pk.txt', 'a')
			#fg = open('l.txt', 'a')
			#f.write(str(pKij_l) + '\n')
			#fg.write(str(l) + '\n')
			#f.close()
			#fg.close()

		self.printPSResults(pKij, larray)
		return pKij,larray
	

	def shearPSCalc(self):

		if self.i > self.shearBins - 1 or self.j > self.shearBins - 1 or self.i < 0 or self.j < 0:
			print "i or j (bin number) is greater than the number of Shear bins: ", self.shearBins
			return

		lc = lensingClassFile.lensingClass(self.h, self.omegamat, self.omegaDE)

		larray = logspace(math.log10(1), math.log10(1e8), self.nummulipoles)
		pKij = []

		F2redshifts = lc.getRedshiftsFromDLS()    	
		#organizedList, numberedList = lc.organizeRedshifts(F2redshifts)
		
		#f_x = interp1d(linspace(0, max(zttiarray), len(Warraytt)), Warraytt)

		#redshifts = linspace(0, max(zttiarray), len(organizedList))
		#Warraytt = f_x(redshifts)
		
		Warrayi = []
		Warrayj = []

		shearbinarray = linspace(0, max(F2redshifts), self.shearBins+1)

		if self.shearBins == 1:

			zttiarray = F2redshifts
			organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
			sortedzttiarray = lc.sort_Redshift_Values(organizedListi)

		    	Warrayi = []
			#print "sortedzttiarray: " + sortedzttiarray
		
		    	for z in sortedzttiarray[0:len(sortedzttiarray) - 1]:
				#print "z                                 : " + str(z)
				Wvali = lc.lensingWeightFunction(z, zttiarray, organizedListi, numberedListi)
				Warrayi.append(Wvali)
		else:

			zttiarray = [ii for ii in F2redshifts if ii >= shearbinarray[self.i] and ii < shearbinarray[self.i+1]]
			zttjarray = [ii for ii in F2redshifts if ii >= shearbinarray[self.j] and ii < shearbinarray[self.j+1]]

			organizedListi, numberedListi = lc.organizeRedshifts(zttiarray)
			sortedzttiarray = lc.sort_Redshift_Values(organizedListi)			
			
			organizedListj, numberedListj = lc.organizeRedshifts(zttjarray)
			sortedzttjarray = lc.sort_Redshift_Values(organizedListj)

		    	Warrayi = []
		    	Warrayj = []
		
		    	for z in sortedzttiarray[1:len(sortedzttiarray) - 1]:
			
				Wvali = lc.lensingWeightFunction(z, sortedzttiarray, organizedListi, numberedListi)
				Wvalj = lc.lensingWeightFunction(z, sortedzttjarray, organizedListj, numberedListj)

				Warrayi.append(Wvali)
				Warrayj.append(Wvalj)
		#print "Size of Warray: ", len(Warrayi), "\nSize of sortedzttiarray: ", len(sortedzttiarray[0:len(Warrayi)])
		lc.plotWZ(sortedzttiarray[0:len(Warrayi)], Warrayi)
			
		cosmoCalcs = cosmocalcs.cosmologyCalculator(lc.h, self.omegamat, self.omegaDE, lc.wX, lc.wXa)
		tf = powerspec.transferFunction(self.omegamat-lc.omegams[2], lc.omegams[2], self.h)
		gro = powerspec.growthFunctionLCDM(self.omegamat, self.omegaDE)
		ps = powerspec.powerSpectrum(tf, gro, self.ns, self.DelRsq)

		ADDM = []
		for z in sortedzttiarray[1:len(zttiarray) - 1]:
			cosmoCalcs.setEmissionRedShift(z)					
			ADDM.append(cosmoCalcs.AngularDiameterDistanceMpc())

		for l in larray:

			pKij_l = 0
			ii = 0

			for z in sortedzttiarray[1:len(zttiarray) - 1]:
				#print "Determining cross-spectra weight function at z: ", z
				#cosmoCalcs.setEmissionRedShift(z)

				D_AMpc = (1 + z)*ADDM[ii]
				Pk = ps.getPk(l/D_AMpc, z) # in (Mpc/h)^3 units

				if self.shearBins == 1:
					pKij_l += Warrayi[ii]*Warrayi[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegamat*pow(1 + z, 3) + (1 - self.omegaDE - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])
				else:
					pKij_l += Warrayi[ii]*Warrayj[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegaDE*pow(1 + z, 3) + (1 - self.omegaDE - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])

				ii += 1

			pKij.append(pKij_l)

			#f = open('Pk.txt', 'a')
			#fg = open('l.txt', 'a')
			#f.write(str(pKij_l) + '\n')
			#fg.write(str(l) + '\n')
			#f.close()
			#fg.close()

		self.printPSResults(pKij, larray)
		return pKij,larray
	
	def printPSResults(self, pKij, larray):
		
		import math
		import matplotlib.pyplot as plt
	
		pKij_l_1 = []
		pKij_l_1_1bins = []
		i = 0
		
		for lt in larray:
			pKij_l_1.append(lt*(lt + 1)*pKij[i])
			pKij_l_1_1bins.append(lt*(lt + 1)*pKij[i]/(2.0*math.pi))
			i += 1
		
		plt.plot(larray, pKij, 'r-')
		plt.semilogy()
		plt.semilogx()
		plt.xlabel('l')
		plt.ylabel('$P_K^{' + str(self.i) + str(self.j) + '}$')
		plt.title('Shear Power Spectrum (' + str(self.shearBins) + ' bins) vs. Multipole Number')
		plt.axis([min(larray), max(larray), min(pKij), max(pKij)])
		plt.grid(True)
		plt.savefig('pKij_l.png')
		plt.close()
		
		plt.plot(larray, pKij_l_1, 'r-')
		plt.semilogy()
		plt.semilogx()
		plt.xlabel('l (multipole number)')
		plt.ylabel('$l(l + 1)/2\pi P_K^{' + str(self.i) + str(self.j) + '}$')
		plt.title('Shear Power Spectrum (' + str(self.shearBins) + ' bins) vs. Multipole Number')
		plt.axis([min(larray), max(larray), min(pKij_l_1), max(pKij_l_1)])
		plt.grid(True)
		plt.savefig('pKij_l_1.png')
		plt.close()
		
class convolutionClass(object):

	def __init__(self, n=1, i=0, j=1, h=0.7, omegamat=0.3, omegaDE=0.7):
		self.i = i
		self.j = j
		self.shearBins = n
		self.h = h
		self.omegamat = omegamat
		self.omegaDE = omegaDE
		self.ns = 0.96
		self.DelRsq = 2.1e-9
		#print 'Initialization of convolutionClass'
	
	def convolutionCalc(self, Pkarray, larray):	

		from scipy.special import jn as J

		bins = 10
		from numpy import linspace		
		limits = linspace(min(larray), max(larray), bins)

		e_plus = []
		e_minus = []

		#theta = linspace(1/60.0*math.pi/180.0, 100.0/60.0*math.pi/180.0, 30)
		theta = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 30)

		P_l = interp1d(larray, Pkarray)

		f_x = lambda l,t: l/(2*math.pi)*J(0, l*t)*P_l(l)
		f_x_m = lambda l,t: l/(2*math.pi)*J(4, l*t)*P_l(l)
		
		for t in theta:

			e_m = 0.0
			e_p = 0.0

			for i in range(bins-1):
				#e_p_, errrop = quad(f_x, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
				e_p_, errrop = quad(f_x, limits[i],  limits[i+1], args=(t,), epsabs = 1e-11, limit = 2000, full_output = False)
				e_p += e_p_

				#e_m_, errrom = quad(f_x_m, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
				e_m_, errrom = quad(f_x_m, limits[i],  limits[i+1], args=(t,), epsabs = 1e-11, limit = 2000, full_output = False)
				e_m += e_m_

			e_plus.append(e_p)
			e_minus.append(e_m)
			

		e_prime = []
		i = 0

		for t in theta:

			e_prime.append(e_minus[i])
			j = 0

			for tt in theta:
	
				if tt > t:

					if j == len(theta) - 1:

						e_prime[i] += 0#4*(e_minus[j]/2.0)*(theta[] - theta[1])/theta[j]
						e_prime[i] -= 0#12 * t**2 * e_minus[j]*(theta[2] - theta[1])/(theta[j]**3)

					else:

						e_prime[i] += 4*(e_minus[j] + e_minus[j+1])/2.0*(theta[j+1] - theta[j])/theta[j]
						e_prime[i] -= 12 * t**2 * (e_minus[j+1] + e_minus[j])/2.0*(theta[j+1] - theta[j])/(theta[j]**3)

				j += 1

			i += 1


		for i in range(len(theta)):

			e_plus[i] = (e_plus[i] + e_prime[i])/2.0
			e_minus[i] = (e_plus[i] - e_prime[i])/2.0

		#thetaDisplay = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 30)
		#e_arrayP = interp1d(e_plus, theta)
		#e_arrayM = interp1d(e_minus, theta)
		#eDisplayplus = e_arrayP(thetaDisplay)
		#eDisplayminus = e_arrayM(thetaDisplay)

		self.printCCresults(e_plus, e_minus, theta)#[180.0*math.pi*60.0*i for i in theta])

		return e_plus, e_minus, theta

	def printCCresults(self, e_plus, e_minus, theta):

		import matplotlib.pyplot as plt

		edat = [0.000174, 0.000239, 0.000146, 0.000157, 0.000145, 0.000112, 0.000094, 0.000124, 0.000121, 0.000129, 0.000112, 0.000111, 0.000093, 0.000082, 0.000058, 0.000043, 0.000026, 0.000015, 0.000010, 0.000009, 0.000005, 0.000006, -0.000007, -0.000012, -0.000018, -0.000018, -0.000013, -0.000011, 0.000004, 0.0]

		eout = open('e_plus.txt', 'w')
		emout = open('e_minus.txt', 'w')

		for i in range(len(theta)):
			theta[i] *= 60.0*180.0/math.pi
			eout.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\n')
			emout.write(str(theta[i]) + '\t' + str(e_minus[i]) + '\n')

		eout.close()
		emout.close()
				
		######################################
		
		plt.plot(theta, e_plus, 'b-')
		plt.plot(theta, edat, 'r*')
		plt.semilogx()
		plt.xlabel('theta')
		plt.ylabel('e_plus/e_minus')
		plt.title('e_plus/e_minus vs. theta')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.axis([1, max(theta), -5e-5, 3e-4])#min(e_plus), max(e_plus)])
		plt.grid(True)
		plt.savefig('e_plus_minus.png')
		plt.close()
		
		"""plt.plot(theta, edat, 'r*')
		plt.semilogx()
		plt.xlabel('theta')
		plt.ylabel('e_minus')
		plt.title('e_minus vs. theta')
		plt.axis([1, max(theta), -5e-5, 3e-4])#min(e_plus), max(e_plus)])
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.grid(True)
		plt.savefig('e_dat.png')
		plt.close()"""
