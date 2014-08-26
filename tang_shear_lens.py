import math
import os.path
import numpy as np
import pyfits

import dls
#import cosmo
#import const



""" Class for computing tangential shear around a lens

Modified from galaxy-galaxy lensing code written by Ami Choi

TODO: write class that uses this class and stacks tangential shear around
      a series of lenses

Last modified: 27 Feb 2014 by AA
"""


class tangShearLens(object):

    """Default tangShearLens object"""
    def __init__(self, ra_lens, dec_lens, z_lens, field_lens, input_srcfile=None, **kws):

        # hard coded parameters
        self.snsq = 0.3*0.3 # shape noise (hard coded)
        
        # read sources
        if input_srcfile is None:
            input_srcfile="Sources.cat"
        self.input_srcfile = input_srcfile
        self._readSources()

        # lens parameters
        self.ra_lens, self.dec_lens = ra_lens, dec_lens
        self.xlens, self.ylens = dls.returnxy(field_lens,ra_lens,dec_lens)
        self.z_lens = z_lens
        #self.distlens = cosmo.D_angular_diameter([self.z_lens], OmegaM, OmegaQ, w0, wa)[0]
        
        # default quality cuts
        demax = 0.25 # maximum ellipticity error
        magmin = 21. # min magnitude
        magmax = 27. # max magnitude
        bmin = 0.4   # min semi-minor axis
        dz = 0.2     # delta-z for cluster foreground cut
        zmax = 1.4   # max photo-z of source
        qcmin = 5    # does not look like a point source (based on 2nd moments)
        usestar = True # use star id to remove star-like objects
        self.set_selection(demax, magmin, magmax, bmin, dz, zmax, qcmin, usestar)
        
        # set shear responsivity correction
        a = 20.
        b = 6e-4
        c = 3.26
        d = 1.036
        self.set_response_correction(a,b,c,d)
        
        # number of bootstrap samples
        self.nboot = 1000
        

    def _readSources(self):
        """ Read in info about the sources """
        hdulist = pyfits.open(self.input_srcfile)
        d = hdulist[1].data
        fileHeader = hdulist[1].header
        ntotal = fileHeader['NAXIS2']
        print ntotal,"sources in file"
        hdulist.close()
        self.xsrc_arr = d.field('x').copy()         # x (in pixels) 
        self.ysrc_arr = d.field('y').copy()         # y (in pixels) 
        self.e1_arr = d.field('e1').copy()          # e1 of all sources 
        self.e2_arr = d.field('e2').copy()          # e2 of all sources 
        self.zp_arr = d.field('z_b').copy()         # photo-z of all sources
        self.de = d.field('de').copy()              # ellipticity error
        self.R_arr = d.field('Rdered').copy()       # R band magnitude
        self.b_arr = d.field('b').copy()            # semi-minor axis
        self.status_arr = d.field('status').copy()  # shear estimation status
        self.flr_arr = d.field('flux_radius').copy()# flux radius
        self.dlsqc_arr = d.field('dlsqc').copy()    # star cut
        self.starid_arr = d.field('starid').copy()  # star id
        self.srcsize = self.xsrc_arr.size           # number of sources in file
        # for testing with compute_shear_other
        #self.rasrc_arr = d.field('RA').copy()
        #self.decsrc_arr = d.field('DEC').copy()


    def _selection(self,minz):
        
        if (self.useStarid):
            good = (self.de<self.demax)&(self.de>0.)&                \
               (self.R_arr>self.magmin)&(self.R_arr<self.magmax)&    \
               (self.zp_arr>=minz)&(self.zp_arr<self.zmax)&          \
               (self.b_arr>self.bmin)&                               \
               (self.status_arr==1)&                                 \
               (self.dlsqc_arr>self.qcmin)&                          \
               (self.starid_arr=='N')        
        else:
            good = (self.de<self.demax)&(self.de>0.)&                \
               (self.R_arr>self.magmin)&(self.R_arr<self.magmax)&    \
               (self.zp_arr<self.zmax)&                              \
               (self.b_arr>self.bmin)&                               \
               (self.status_arr==1)&                                 \
               (self.zp_arr>=minz)&                                  \
               (self.dlsqc_arr>self.qcmin)                               
   
        return good


    def _response_correction(self):
        """ Shear responsivity correction """
        
        x = (self.R_arr - self.a)
        x[x<0.] = 0. # set all negative values to zero
        mgmma = self.b*pow(x, self.c) + self.d
        return mgmma


    def set_selection(self, demax, magmin, magmax, bmin, dz, zmax, qcmin, usestar):
        """ If want to change these selection parameters """
        
        self.demax = demax       # maximum ellipticity error
        self.magmin = magmin     # min magnitude
        self.magmax = magmax     # max magnitude
        self.bmin = bmin         # min semi-minor axis
        self.dz = dz             # foreground cut
        self.zmax = zmax         # max photo-z
        self.qcmin = qcmin       # min qc
        self.useStarid = usestar # use star id field


    def set_response_correction(self,a,b,c,d):
        """ Set responsivity correction parameters, to switch off correction
            set_response_correction(0.,1.,0.,0.)                             """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        

    def compute_shear(self,ri,rf):
        """ Compute tangential shears and store in big array
        ri and rf are the min and max radii in arcsec 
        dz is the photo-z foreground cut """
        
        # source-lens separation in pixels
        delta = np.hypot(self.xsrc_arr-self.xlens,self.ysrc_arr-self.ylens)
        # now in arcsecs
        delta *= dls.pixsc

        # find all sources with ri<= separation <=rf and with z>z_lens+deltaz 
        # and that match the other quality cuts
        minz = self.z_lens + self.dz
        good = (delta<=rf)&(delta>=ri)&self._selection(minz)
        e1_src = self.e1_arr[good].copy()
        e2_src = self.e2_arr[good].copy() # why was it *-1? cfhtls coord flip!
        dy_arcsec = (self.ysrc_arr[good]-self.ylens)*dls.pixsc # arcsec
        dx_arcsec = (self.xsrc_arr[good]-self.xlens)*dls.pixsc # arcsec
        theta = delta[good]
        
        # calculate cos2phi = (dx/dr)^2 - (dy/dr)^2
        cosarray = (dx_arcsec*dx_arcsec - dy_arcsec*dy_arcsec)
        cosarray = cosarray/(theta*theta)
        # calculate sin2phi = 2dydx/dr^2
        sinarray = (2.*dx_arcsec*dy_arcsec)/(theta*theta)
        # tangential and cross shear
        mgmma = self._response_correction()
        #print "min,max mgamma =",np.amin(mgmma[good]),np.amax(mgmma[good])
        e_t = mgmma[good]*(-e1_src*cosarray-e2_src*sinarray)
        e_x = mgmma[good]*(+e1_src*sinarray-e2_src*cosarray)
        # weight
        wt = 1./(self.de[good]*self.de[good] + self.snsq)
        
        
        return theta, e_t, e_x, wt, self.zp_arr[good]


    def linearbin_shear(self, theta, e_t, e_x, wt, bmin, bmax, nbin, weight=True):
        """ Radially bin the tangential shear with linear spaced bins

            theta is source-lens separation in arcsecs
            e_t, e_x are tangential and cross shear components
            wt is the weight of each shear (1/sqrt(sig_e^2+sig_sn^2)
            bmin, bmax are the bin min and max in arcsecs
            nbin is the number of radial bins   

            need to implement shear weigting          
        """

        number_of_sources = len(theta)
        binsize = (bmax-bmin)/nbin

        sum_rad = np.zeros([nbin])
        sum_etd = np.zeros([nbin])
        sum_enull = np.zeros([nbin])
        sumsq_etd = np.zeros([nbin])
        sumsq_enull = np.zeros([nbin])
        n_in_bin = np.zeros([nbin])
        for i in range(number_of_sources):

            idgal = math.floor((theta[i] - bmin)/binsize)

            if (idgal>=0 and idgal<nbin):
                sum_rad[idgal] += theta[i]
                sum_etd[idgal] += e_t[i]
                sum_enull[idgal] += e_x[i]
                sumsq_etd[idgal] += (e_t[i]*e_t[i])
                sumsq_enull[idgal] += (e_x[i]*e_x[i])
                n_in_bin[idgal] += 1

        print n_in_bin.sum(),"galaxies out of",number_of_sources,"made it into bins"
        mean_rad = sum_rad/n_in_bin
        mean_et = sum_etd/n_in_bin
        mean_ex = sum_enull/n_in_bin
        # this is not what the error should be below, should bootstrap
        err_et = np.sqrt((sumsq_etd - mean_et)/n_in_bin/n_in_bin)
        err_ex =  np.sqrt((sumsq_enull - mean_ex)/n_in_bin/n_in_bin)
        return mean_rad, mean_et, mean_ex, err_et, err_ex, n_in_bin


    def logbin_shear(self, theta, e_t, e_x, wt, bmin, bmax, nbin, weight=True):
        """ Radially bin the tangential shear with log spaced bins
            
            theta is source-lens separation in arcsecs
            e_t, e_x are tangential and cross shear components
            wt is the weight of each shear (1/sqrt(sig_e^2+sig_sn^2)
            bmin, bmax are the bin min and max in arcsecs
            nbin is the number of radial bins   
        
        """
      
        number_of_sources = len(theta)
        logbinsize = (math.log(bmax)-math.log(bmin))/nbin
        
        sum_rad = np.zeros([nbin])
        sum_etd = np.zeros([nbin])
        sum_enull = np.zeros([nbin])
        sumsq_etd = np.zeros([nbin])
        sumsq_enull = np.zeros([nbin])
        sum_wt = np.zeros([nbin])
        n_in_bin = np.zeros([nbin])
        
        for i in range(number_of_sources):
            #print 'chk: i =',i,'theta[i] =',theta[i],'log(bmin) =',math.log(bmin),
            #print 'logbinsize =',logbinsize
            
            idgal = math.floor((math.log(theta[i]) - math.log(bmin))/logbinsize)
            
            if (idgal>=0 and idgal<nbin):
                sum_rad[idgal] += theta[i]
                if weight:
                    sum_etd[idgal] += (e_t[i]*wt[i])
                    sum_enull[idgal] += (e_x[i]*wt[i])
                else:
                    sum_etd[idgal] += e_t[i]
                    sum_enull[idgal] += e_x[i]
                sumsq_etd[idgal] += (e_t[i]*e_t[i])
                sumsq_enull[idgal] += (e_x[i]*e_x[i])
                n_in_bin[idgal] += 1
                sum_wt[idgal] += wt[i]
        
        print int(n_in_bin.sum()),"galaxies out of",number_of_sources,"made it into bins"
        mean_rad = sum_rad/n_in_bin
        
        if weight:
            print "Performing weighted mean in bins"
            mean_et = sum_etd/sum_wt
            mean_ex = sum_enull/sum_wt
            err_et = np.sqrt(1./sum_wt)
            err_ex =  np.sqrt(1./sum_wt)
        else:
            print "Performing unweighted mean in bins"
            mean_et = sum_etd/n_in_bin
            mean_ex = sum_enull/n_in_bin
            # this is not what the error should be below, should bootstrap
            err_et = np.sqrt((sumsq_etd - mean_et)/n_in_bin/n_in_bin)
            err_ex =  np.sqrt((sumsq_enull - mean_ex)/n_in_bin/n_in_bin)
            
        return mean_rad, mean_et, mean_ex, err_et, err_ex, n_in_bin


    def _write_header(self, f, nsources, write_shear_prof=True):
        """ Write header containing info of how sources were selected for 
            tangential shear profile
        """
        if (write_shear_prof):
            header = "# File sources read from " + str(self.input_srcfile)+"\n"
            f.write(header)
            header = "# Cluster ra,dec,z = " + str(self.ra_lens) + ","+str(self.dec_lens)
            header +="," + str(self.z_lens) + "\n"
            f.write(header)
        header = "# number of sources = " + str(int(nsources)) + "\n"
        f.write(header)
        header = "# de < " + str(self.demax) + " (maximum ellipticity error)\n"
        f.write(header)
        header = "# " + str(self.magmin) + " < mag < " + str(self.magmax) 
        header +=" (magnitude range) \n"
        f.write(header)
        header = "# b > " + str(self.bmin) + " (min semi-minor axis in pixels)\n"
        f.write(header)
        header = "# " + str(self.z_lens+self.dz) + " < zphot < " + str(self.zmax)
        header +=" (photo-z range of sources)\n"
        f.write(header)        header = "# qc > " + str(self.qcmin)
        if (self.useStarid):
            header += ", and starid = null \n"
        else:
            header += "\n"
        f.write(header)


    def write_shear_profile(self, filename, rad, tang, cross, errt, errx, nb):
        """ write the tangential shear profile tang(rad) and cross(rad) to a 
            file """
        
        # check if file exists
        if (os.path.isfile(filename)):
            print "ERROR! file ", filename ," exists!"
            return
            
        # write file
        nsources = nb.sum()
        nbin = len(rad)
        f = open(filename,'w')
        self._write_header(f, nsources)
        
        header = "# Columns: radius (arcsecs) - gamma_t - gamma_x - e_gamma_t"
        header +=" - e_gamma_x - n_in_bin \n"
        f.write(header)
        for i in range(nbin):
            line = str(rad[i]) + "  " + str(tang[i]) + "  " +  str(cross[i]) \
            + "  " + str(errt[i]) + "  " +  str(errx[i]) + "  " + str(int(nb[i])) +"\n"
            f.write(line)
        f.close()
        
        
    def write_sigma_profile(self, filename, rad, dsigma, err, nb):
        """ write the excess surface mass density profile to a file """
        
        # check if file exists
        if (os.path.isfile(filename)):
            print "ERROR! file ", filename ," exists!"
            return
            
        # write file
        nsources = 0
        for n in nb:
            nsources+=n
        nbin = len(rad)
        f = open(filename,'w')
        self._write_header(f, nsources, False)
        
        header = "# Columns: radius - delta_Sigma - error n_in_bin \n"
        f.write(header)
        for i in range(nbin):
            line = str(rad[i]) + "  " + str(dsigma[i]) + "  " +  str(err[i]) \
            + "  " + str(int(nb[i])) +"\n"
            f.write(line)
        f.close()


    def set_nboot(nboot):
        self.nboot = nboot
        
        
# Not part of tangShearLens 

def stack_dsigma(rad, dSigma, weight, bmin, bmax, nbin, nboot=1000):
    """ Stack the excess surface mass density in log spaced bins and 
        calculate the bootstrapped error"""
        
    # bootstrap average and error in bin
    out_t = []
    out_t_boot = []      
    # Calculate the radial bin boundaries
    bin_centers = []
    bin_edges = np.zeros((nbin+1))
    for i in range(nbin+1):
        bin_edges[i] = 10**(np.log10(bmin)+(i*(np.log10(bmax/bmin))/nbin))
        
    # number of sources in each bins
    n_in_bins = []
    
    print 'min rad =',min(rad),'max rad =',max(rad)
    print weight
        
    for i in range(nbin):
        
        print "Bin %1d/%1d" %(i+1,nbin),':',
        
        # find all sources in current bin
        bin_low,bin_up = bin_edges[i],bin_edges[i+1]
        print bin_low,'to',bin_up
        cond = ((rad<bin_up)&(rad>=bin_low))
        #print cond[:10]
        dSigma_bin = dSigma[cond].copy()
        print len(dSigma_bin),'in bin'
        weight_bin = weight[cond].copy()
        

        # initialize bootstrap sample to zero
        boot_dsbin = np.empty(0)

        for j in range(nboot):
            
            # randomly select data points from data in bin
            # (could have repetitions)
            choices = np.random.random_integers(0,dSigma_bin.size-1,dSigma_bin.size)
            boot_ds = dSigma_bin[choices].copy()
            boot_weight = weight_bin[choices].copy()
                
            # append average to previous bootstrap samples of this bin
            boot_dsbin = np.append(boot_dsbin,np.average(
                                   boot_ds, weights=boot_weight))

        # weighted average of dSigma in the bin:
        out_t.append(np.average(dSigma_bin,weights=weight_bin))
        # bootstrap error of dSigma in bin
        out_t_boot.append(boot_dsbin.std())

        # central value of bin
        bin_centers.append((bin_low+bin_up)/2.)
            
        # add to total in each bin
        n_in_bins.append(dSigma_bin.size)
            
    return bin_centers, out_t, out_t_boot, n_in_bins

