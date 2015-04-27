import sys
import lensingClassFile
import cosmocalcs
from numpy import linspace
from numpy import logspace
from scipy.interpolate import interp1d
from scipy.integrate import quad as sciquad
import mpmath
from mpmath import quad
from mpmath import quadosc
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
		self.nummulipoles = 5000
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

		larray = logspace(math.log10(1), math.log10(1e6), self.nummulipoles)
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
					#pKij_l += Warrayi[ii]*Warrayi[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegamat*pow(1 + z, 3) + (1 - self.omegaDE - lc.omegams[2])*pow(1 + z, 2) + lc.omegams[2], 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])
					pKij_l += Warrayi[ii]*Warrayi[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegamat*pow(1 + z, 3) + (1 - self.omegaDE - self.omegamat)*pow(1 + z, 2) + self.omegaDE, 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])
				else:
					pKij_l += Warrayi[ii]*Warrayj[ii]*Pk/(D_AMpc*D_AMpc*100*self.h*pow(self.omegaDE*pow(1 + z, 3) + (1 - self.omegaDE - self.omegamat)*pow(1 + z, 2) + self.omegaDE, 0.5))*(sortedzttiarray[1] - sortedzttiarray[0])

				ii += 1

			pKij.append(pKij_l)

			#f = open('Pk.txt', 'a')
			#fg = open('l.txt', 'a')
			#f.write(str(pKij_l) + '\n')
			#fg.write(str(l) + '\n')
			#f.close()
			#fg.close()

		#self.printPSResults(pKij, larray)
		return pKij,larray
	

	def shearPSCalc(self):

		if self.i > self.shearBins - 1 or self.j > self.shearBins - 1 or self.i < 0 or self.j < 0:
			print "i or j (bin number) is greater than the number of Shear bins: ", self.shearBins
			return

		lc = lensingClassFile.lensingClass(self.h, self.omegamat, self.omegaDE)

		larray = logspace(math.log10(1), math.log10(1e6), self.nummulipoles)
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

		#self.printPSResults(pKij, larray)
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

	def taylorExp(self, x):
		return x - x**3/18 + x**5/600 - x**7/35280 + x**9/3265920 - x**11/439084800 

	def printBessel(self, Plarray, larray):

		from scipy.special import jn as J
		from numpy import logspace
		theta = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 30)
		import matplotlib.pyplot as plt

		t = theta[0]
		y = []
		for l in larray:
			y.append(J(0,t*l))	

		plt.plot(larray, y, 'r-')
		plt.xlabel('l')
		plt.ylabel('J_0(l*t)')
		plt.title('J_0 vs. l')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.axis([1, max(larray)/100, .1, -.1])#min(e_plus), max(e_plus)])
		plt.grid(True)

		t = theta[len(theta)-1]
		y = []
		for l in larray:
			y.append(J(0,t*l))	

		plt.plot(larray, y, 'b-')
		plt.xlabel('l')
		plt.ylabel('J_0(l*t)')
		plt.title('J_0 vs. l')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.axis([1, max(larray)/100, .1, -.1])#min(e_plus), max(e_plus)])
		plt.grid(True)
			
		plt.savefig('J_l.png')
		plt.close()



		t = theta[0]
		y = []
		for l in larray:
			y.append(J(0,t*l))	

		plt.plot(larray, y, 'r-')
		plt.xlabel('l')
		plt.ylabel('J_0(l*t)')
		plt.title('J_0 vs. l')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.axis([1, 10000, 1, -1])#min(e_plus), max(e_plus)])
		plt.grid(True)

		t = theta[len(theta)-1]
		y = []
		for l in larray:
			y.append(J(0,t*l))	

		plt.plot(larray, y, 'b-')
		plt.xlabel('l')
		plt.ylabel('J_0(l*t)')
		plt.title('J_0 vs. l')
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.axis([1, 10000, 1, -1])#min(e_plus), max(e_plus)])
		plt.grid(True)
			
		plt.savefig('J_l2.png')
		plt.close()




	def convolutionCalc(self, Pkarray, larray):	

		print 'Convolution begineth'

		from scipy.special import jn as J
		#f = lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)
		#g = lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x)


		bins = 30
		from numpy import linspace

		limits = linspace(min(larray), max(larray), bins)
		

		e_plus = []
		e_minus = []

		theta = logspace(math.log10(1.0/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 30)
		#theta = logspace(math.log10(0.01/60.0*math.pi/180.0), math.log10(100.0/60.0*math.pi/180.0), 3000)

		P_l = interp1d(larray, Pkarray)

		#f = open('P_l.txt', 'w')
		#for l in larray:
		#	f.write(str(l) + '\t' + str(P_l(l)) + '\n')
		#f.close()


###########################################################################################
# Analytical attempt
#
#		f_x = lambda l,t: l/(2*math.pi)*J(0, l*t)*P_l(l)
#		f_x_m = lambda l,t: l/(2*math.pi)*J(4, l*t)*P_l(l)		
#
#		for t in theta:
#
#			i = 0
#			e_p = 0.0
#			e_m = 0.0
#
#			e_p_ = 0.0
#			e_m_ = 0.0
#
#			while (i+1)*3.1415926535/t <= max(larray):	
#				if (i+1)*3.1415926535/t >= 12000.0 and i*3.1415926535/t >= 1.0:
#					#print 'P_l', i*3.1415926535/t
#					df = (P_l((i+1)*3.1415926535/t) - P_l(i*3.1415926535/t))/((i+1)*3.1415926535/t - i*3.1415926535/t) 		
#					e_p_ = (P_l(i*math.pi/t) - df*i*math.pi/t) * (self.taylorExp((i+1)*math.pi) - self.taylorExp(i*math.pi)) - df/t * (math.cos((i+1)*math.pi) - math.cos(i*math.pi))
#					
#				else:
#					print 'tehta', 20000.0*t/3.1415926535
#
#					e_p_, errrop = quad(f_x, 1,  12000.0*t/3.1415926535, args=(t,), epsabs = 1e-10, limit = 10000, full_output = False)
#					e_m_, errrom = quad(f_x_m, 1,  12000.0*t/3.1415926535, args=(t,), epsabs = 1e-10, limit = 10000, full_output = False)			
#				
#				e_p += 1/(2*math.pi)*e_p_
#				e_m += 1/(2*math.pi)*e_m_
#				i += 1
#
#
#			e_plus.append(e_p)
#			#e_minus.append(e_m)
#
###########################################################################################

###########################################################################################
#  Quadosc
#
#		f = lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x)
#		g = lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x)
#
#		for t in theta:
#
#			#i = -1
#			#while lt <= 1e6:	
#			#	i += 1
#			#	limits.append(i*3.1315926535/t)
#			#	lt = limits[i]
#			#	print lt
#
#			print 'Theta', t*180.0/3.14*60
#				
#
#			e_m = 0.0
#			e_p = 0.0
#
#			for i in range(bins-1):
#
#				e_p_ = mpmath.quadosc(f, [limits[i], limits[i+1]], omega=t, verbose=True, maxdegree=20, method='gauss-legendre')				
#				e_p += e_p_
#
#				e_m_ = mpmath.quadosc(g, [limits[i], limits[i+1]], omega=t, verbose=True, maxdegree=20, method='gauss-legendre')				
#				e_m += e_m_
#
#			e_plus.append(e_p)
#			e_minus.append(e_m)
#
###########################################################################################

###########################################################################################
# quad method
#
#		f1 = open('testingshit.txt', 'w')
#		f1.close()
#
#		#self.printBessel(Pkarray, larray)
#
#		f = lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x)
#		g = lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x)
#
#
#		for t in theta:
#
#			#i = -1
#			#while lt <= 1e6:	
#			#	i += 1
#			#	limits.append(i*3.1315926535/t)
#			#	lt = limits[i]
#			#	print lt
#
#			print 'Theta', t*180.0/3.14*60
#				
#
#			e_m = 0.0
#			e_p = 0.0
#
#			for i in range(bins-1):
#				#e_p_, errrop = quad(f_x, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
#				#e_p_, errrop = quad(f_x, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 10000, full_output = False)
#				e_p_ = mpmath.quad(lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x), [limits[i], limits[i+1]], maxdegree=4)#, method='gauss-legendre')#, verbose=True)
#				e_p += float(e_p_)
#
#				#e_m_, errrom = quad(f_x_m, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
#				#e_m_, errrom = quad(f_x_m, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 20000, full_output = False)
#				e_m_ = mpmath.quad(lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x), [limits[i], limits[i+1]], maxdegree=4)#, method='gauss-legendre')#, verbose=True)				
#				e_m += float(e_m_)
#
#			e_plus.append(e_p)
#			e_minus.append(e_m)
#
#		i = 0
#		f1 = open('testingshit.txt', 'a')
#		for t in theta:
#			f1.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\t' + str(e_minus[i]) + '\n')			
#			i += 1
#		f1.close()
#
#
###########################################################################################

###########################################################################################
# original method

		f_x = lambda l,t: l/(2*math.pi)*J(0, l*t)*P_l(l)
		f_x_m = lambda l,t: l/(2*math.pi)*J(4, l*t)*P_l(l)

		for t in theta:

			#i = -1
			#while lt <= 1e6:	
			#	i += 1
			#	limits.append(i*3.1315926535/t)
			#	lt = limits[i]
			#	print lt

			print 'Theta', t*180.0/3.14*60
				

			e_m = 0.0
			e_p = 0.0

			for i in range(bins-1):
				#e_p_, errrop = quad(f_x, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
				e_p_, errrop = sciquad(f_x, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 10000, full_output = False)
				#e_p_ = mpmath.quad(lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x), [limits[i], limits[i+1]], maxdegree=2)#, method='gauss-legendre')#, verbose=True)
				e_p += e_p_

				#e_m_, errrom = quad(f_x_m, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
				e_m_, errrom = sciquad(f_x_m, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 10000, full_output = False)
				#e_m_ = mpmath.quad(lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x), [limits[i], limits[i+1]], maxdegree=2)#, method='gauss-legendre')#, verbose=True)				
				e_m += e_m_

			e_plus.append(e_p)
			e_minus.append(e_m)

		#i = 0
		#f1 = open('testingshit.txt', 'a')
		#for t in theta:
		#	f1.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\t' + str(e_minus[i]) + '\n')			
		#	i += 1
		#f1.close()

################################################################################################

################################################################################################
# First and then every nth integration
#
#
#		f = lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x)
#		g = lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x)
#
#		for t in theta:
#
#			#i = -1
#			#while lt <= 1e6:	
#			#	i += 1
#			#	limits.append(i*3.1315926535/t)
#			#	lt = limits[i]
#			#	print lt
#
#			print 'Theta', t*180.0/3.14*60
#				
#			e_m = 0.0
#			e_p = 0.0
#
#			e_p_, errrop = quad(f, min(larray),  10800, args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
#			e_p += e_p_		
#			e_m_, errrom = quad(f_x_m, min(larray),  5*10800, args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)				
#			e_m += e_m_
#
#			
#
#
#			for i in range(bins-1):
#				#e_p_, errrop = quad(f_x, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
#				#e_p_, errrop = quad(f_x, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 10000, full_output = False)
#				e_p_ = mpmath.quad(lambda x: x/(2*math.pi)*math.sin(x*t)/(x*t)*P_l(x), [limits[i], limits[i+1]], verbose=True, maxdegree=100)
#				e_p += e_p_
#
#				#e_m_, errrom = quad(f_x_m, min(larray),  max(larray), args=(t,), epsabs = 1e-10, limit = 3000, full_output = False)
#				#e_m_, errrom = quad(f_x_m, limits[i],  limits[i+1], args=(t,), epsabs = 1e-9, limit = 20000, full_output = False)
#				e_m_ = mpmath.quad(lambda x: x/(2*math.pi)*((105/(x*t)**4 - 45/(x*t)**2 + 1)*math.sin(x*t)/(x*t) + (-105/(x*t)**2 + 10)*math.cos(x*t)/(x*t)**2)*P_l(x), [limits[i], limits[i+1]], verbose=True, maxdegree=100)				
#				e_m += e_m_
#
#			e_plus.append(e_p)
#			e_minus.append(e_m)
#
################################################################################################

################################################################################################
# Integration by period
#
#		for t in theta:
#
#			e_m = 0.0
#			e_p = 0.0
#			i = -1
#
#			l = 0
#
#			print 'Theta', 180.0/3.1415926535*t
#
#			while (i+2)*3.1415926535/t <= max(larray):	
#
#				i += 1
#
#				#print 'There', (i)*3.1415926535/t
#
#				if (i)*3.1415926535/t >= 1.0:
#
#					#print 'Here', (i)*3.1415926535/t
#
#					#f.write(str(l) + '\t' + str(t) + '\n')
#
#				
#					#l = (i+1)*3.1315926535/t
#		
#					e_p_, errrop = quad(f_x, i*3.1415926535/t,  (i+1)*3.1415926535/t, args=(t,), epsabs = 1e-10, limit = 2, full_output = False)
#					e_p += e_p_
#
#					e_m_, errrom = quad(f_x_m, i*3.1415926535/t,  (i+1)*3.1415926535/t, args=(t,), epsabs = 1e-10, limit = 1, full_output = False)
#					e_m += e_m_
#
#			e_plus.append(e_p)
#			e_minus.append(e_m)
#
#################################################################################################
			

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


		#f.close()

		for i in range(len(theta)):

			e_minus[i] = (e_plus[i] - e_prime[i])/2.0
			e_plus[i] = (e_plus[i] + e_prime[i])/2.0


		self.printCCresults(e_plus, e_minus, theta)#[180.0*math.pi*60.0*i for i in theta])

		#print "HERE"
		#for t in theta:
		#	print t, e_plus, e_minus

		return e_plus, e_minus, theta

	def printCCresults(self, e_plus, e_minus, theta):

		import matplotlib.pyplot as plt

		edat = [0.000174, 0.000239, 0.000146, 0.000157, 0.000145, 0.000112, 0.000094, 0.000124, 0.000121, 0.000129, 0.000112, 0.000111, 0.000093, 0.000082, 0.000058, 0.000043, 0.000026, 0.000015, 0.000010, 0.000009, 0.000005, 0.000006, -0.000007, -0.000012, -0.000018, -0.000018, -0.000013, -0.000011, 0.000004, 0.0]

		e_low = [5.54E-05, 2.88E-05, 1.63E-05, 2.81E-06, 1.33E-05, 1.84E-05, 2.13E-05, 4.72E-06, 9.09E-06, 1.55E-05, 1.27E-05, 5.85E-06, 8.85E-06, 1.26E-05, 5.29E-06, 7.18E-06, 4.48E-06, 7.37E-06, 5.60E-06, 5.71E-06, 7.18E-06, 3.13E-06, 4.31E-06, 2.67E-06, 2.12E-06, 5.40E-06, 4.15E-06, 2.51E-06, 1.60E-06, 0.0]

		e_high = [4.60E-06, 3.40E-05, 5.25E-05, 4.22E-05, 9.37E-06, 3.13E-05, 1.84E-05, 1.77E-05, 1.46E-05, 1.27E-05, 1.03E-05, 1.63E-05, 7.00E-06, 1.17E-05, 1.32E-05, 1.14E-05, 9.26E-06, 8.24E-06, 1.08E-05, 5.05E-06, 6.98E-06, 2.24E-06, 2.15E-06, 1.20E-06, 2.80E-06, 2.86E-06, 3.41E-06, 3.56E-06, 7.33E-07, 0.0]


		#eout = open('e_plus.txt', 'w')
		#emout = open('e_minus.txt', 'w')

		for i in range(len(theta)):
			theta[i] *= 60.0*180.0/math.pi
		#	eout.write(str(theta[i]) + '\t' + str(e_plus[i]) + '\n')
		#	emout.write(str(theta[i]) + '\t' + str(e_minus[i]) + '\n')
		#
		#eout.close()
		#emout.close()
				
		######################################
		
		plt.plot(theta, e_plus, 'r-')
		#plt.plot(theta, edat, 'ro')
		plt.errorbar(theta, edat, yerr=[e_low, e_high], fmt='-o')		
		plt.semilogx()
		plt.xlabel('theta (arcmin)')
		plt.ylabel('e_E')
		plt.title('e_E vs. theta')
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
