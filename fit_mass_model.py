# 
import os.path
import math
import cmath
import const
import cosmocalcs
import shear
import numpy as np
import time

def interp(xvals,yvals,x):
    
    dx = xvals[1] - xvals[0]
    xmin = xvals[0]
    imax = len(xvals)-1
    
    k = long(math.floor((x-xmin)/dx))
    
    if (k < 1):
        y = yvals[0] + (x-xvals[0])*(yvals[1]-yvals[0])/(xvals[1]-xvals[0])
    elif (k > imax):
        y = yvals[imax] + (x-xvals[imax])*(yvals[imax]-yvals[imax-1])/(xvals[imax]-xvals[imax-1])
    else:
        y = yvals[k-1] + (x-xvals[k-1])*(yvals[k]-yvals[k-1])/(xvals[k]-xvals[k-1])
    
    return y


class fitMassModel(object):

    """Default fitNFWModel object"""
    def __init__(self, z_lens, data, cc):
        
        # NFW model
        self.c = 5
        self.M200min = 0.1  # in units of 1e14
        self.M200max = 10.
        self.nM200 = 35
        
        # SIS model
        self.sigmamin = 0. # in units of km/s (velocity dispersion)
        self.sigmamax = 1500.
        self.nsigma = 60
        
        # cluster
        self.z_lens = z_lens
        
        # data, a numpy array with columns: theta, gt, et, z
        self.data = data
        
        # cosmology
        #self.cosmology = {}
        #self.cosmology["OmegaM"] = 0.3
        self.cosmology = cc
        
        
    def set_concentration(self,c):
        self.c = c
        
        
    def set_cosmology(self,cc):
        self.cosmology = cc
        
        
    def set_m200grid(self,M200min, M200max, nM200):
        self.M200min = M200min
        self.M200max = M200max
        self.nM200 = nM200
        
        
    def set_sigmagrid(self,sigmamin, sigmamax, nsigma):
        self.sigmamin = sigmamin
        self.sigmamax = sigmamax
        self.nsigma = nsigma
        
        
        
    def fitNFW(self):
        dM200 = (self.M200max - self.M200min)/(self.nM200-1) 
        
        
        # set chisq array to be size nM200
        chisq_array = np.zeros([self.nM200])
        
        for i in range(self.nM200):
            print "Fitting mass",i+1,"of",self.nM200
            M200 = 1e14*(self.M200min + dM200*i)
            nfwpars = {}
            nfwpars["M200"] = M200
            nfwpars["c"] = self.c
        
            chisq = 0
            for datum in self.data:
                
                theta_g = datum[0]
                gt_g = datum[1]
                err_g = datum[2]
                z_g = datum[3]
                
                #start_time = time.time()
                gamma_p,kappa_p = shear.nfwshear(theta_g, nfwpars, self.z_lens, z_g, self.cosmology)
                #end_time = time.time()
                #print "time to calculate one predicted shear =",end_time-start_time # =0.2s 
                #print gamma_p
                gt_p = gamma_p/(1.-kappa_p)
                
                chisq += ((gt_p - gt_g)*(gt_p - gt_g))/(2*err_g*err_g)
                
            print chisq
            chisq_array[i] = chisq
            #print chisq_array[i]
            
        return chisq_array
        
        
    def fitSIS(self):
        dsigma = (self.sigmamax - self.sigmamin)/(self.nsigma-1) 

        # set chisq array to be size nM200
        chisq_array = np.zeros([self.nsigma])
        
        for i in range(self.nsigma):
            print "Fitting mass",i+1,"of",self.nsigma
            sigma = self.sigmamin + dsigma*i
        
            chisq = 0
            for datum in self.data:
                
                theta_g = datum[0]
                gt_g = datum[1]
                err_g = datum[2]
                z_g = datum[3]
                
                start_time = time.time()
                gamma_p,kappa_p = shear.sisshear(theta_g, sigma, self.z_lens, z_g, self.cosmology)
                end_time = time.time()
                print "time to calculate one predicted shear =",end_time-start_time # =0.2s 
                print gamma_p
                gt_p = gamma_p #/(1.-kappa_p)
                
                chisq += ((gt_p - gt_g)*(gt_p - gt_g))/(2*err_g*err_g)
                
            print chisq
            chisq_array[i] = chisq
            #print chisq_array[i]
            
        return chisq_array
        
        
    def write_chisq(self,outfile,chisq_array):
    
        # check if file exists
        if (os.path.isfile(outfile)):
            print "ERROR! file ", outfile ," exists!"
            return
    
        f = open(outfile,'w')

        dM200 = (self.M200max - self.M200min)/(self.nM200-1) 
        for i in range(self.nM200):
            M200 = 1e14*(self.M200min + dM200*i)
            
            line = str(M200[i]) + "  " + str(chisq_array[i]) + "\n"
            f.write(line)
        f.close()
        
        
