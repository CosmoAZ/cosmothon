## @package simclusterfield
# 
# this module really for simulating clusters
# will want to rename
#
# [More details.]

import math
import os.path
import numpy as np
import scipy.integrate as integ
import random

import const
import cosmocalcs
import galdist


        
## Class for computing cluster number density profile
#
# Last modified: 9 July 2014 by AA
class plummerProfile(object):

    ## The constructor
    # @param zcluster    redshift of cluster
    # @param rcore       cluster core radius (where profile changes) in Mpc
    # @param nrcut       radius of cluster edge is at nrcut*rcore (where profile set to zero)
    def __init__(self, zcluster, rcore=0.150, nrcut=10.):
        
        self.zcluster = zcluster
        self.rcore = rcore  # Mpc
        self.rcut = nrcut*self.rcore  # Mpc
        self.ccalc = cosmocalcs.cosmologyCalculator()
        self.ccalc.setEmissionRedShift(self.zcluster)
        
            
    ## Return Plummer profile
    # @param r    radius in Mpc
    def getProfile(self, r):
            
        P=0.
        if (r<self.rcut):
            P = pow(1. + (r/self.rcore)*(r/self.rcore),-0.5) - \
                pow(1. + (self.rcut/self.rcore)*(self.rcut/self.rcore),-0.5)
        return P


    ## Return Plummer profile
    # @param theta    angular radius in radians
    def getProfileRadians(self, theta):
            
        da = self.ccalc.AngularDiameterDistanceMpc()
        r = da*theta
        P = getProfile(r)
        return P


## Class for drawing object redshift and position
class simObject(object):

    ## The constructor
    # @param ramin     minimum ra to draw
    # @param ramax     maximum ra to draw
    # @param decmin    minimum dec to draw
    # @param decmax    maximum dec to draw
    # @param zlow      minimum redshift to draw
    # @param zhigh     maximum redshift to draw
    # @param nofz      analytic n(z) function
    def __init__(self,ramin, ramax, decmin, decmax, zlow, zhigh, nofz):
    
        self.ramin = ramin
        self.ramax = ramax
        self.decmin = decmin
        self.decmax = decmax
        self.zlow = zlow
        self.zhigh = zhigh
        self.nofz = nofz
        self.zAtMax = self.nofz.zAtMaxNofz()
    
    
    ## Draw object redshift
    def drawRedshift(self):
    
        zobj = -1.
        nrand = -1. 
        while (True):
            zobj = random.uniform(self.zlow,self.zhigh)
            nrand = random.uniform(0., self.zAtMax)

            nz = self.nofz.getNofZ(zobj)
            if (nrand <= nz):
                break
                
        return zobj
        

    ## Draw a object redshift and angular position
    # returns: redshift, ra, dec
    def drawObject(self):

        # draw galaxy redshift
        zobj = self.drawRedshift()

        # draw galaxy angular position
        raobj = random.uniform(self.ramin, self.ramax)
        decobj = random.uniform(self.decmin, self.decmax)
       
        return zobj,raobj,decobj
        
        
    ## Print sim info
    def printInfo(self):

        print 'Object properties:'
        print 'Simulating objects over right ascensions',self.ramin,'to',self.ramax,
        print ', declinations',self.decmin,'to',self.decmax,'and from redshift',
        print self.zlow,'to',self.zhigh
    
        
## Class to simulate unclustered background galaxies
#
# Last modified 21 July 2014
class simBackgroundGals(simObject):

    ## The constructor
    # @param nofz          analytic n(z) function
    # @param neff          observed number of galaxies per square arcminute
    # @param field_center  position in ra, dec (degrees) of field center
    # @param field_size    area of field in square degrees (assumed to be square)
    # @param zrange        simulate galaxies within redshift range zrange=zlow,zhigh
    def __init__(self, nofz, neff=50, field_center=(135.,30.), field_size=4., zrange=(1.,2.)):
         
        # total number of galaxies in the field
        self.ntotal = int(math.floor(neff*field_size*3600.))  

        # init class that does the object position simulation
        ramin = field_center[0] - math.sqrt(field_size)/2.
        ramax = field_center[0] + math.sqrt(field_size)/2.
        decmin = field_center[1] - math.sqrt(field_size)/2.
        decmax = field_center[1] + math.sqrt(field_size)/2.
        simObject.__init__(self, ramin, ramax, decmin, decmax, zrange[0], zrange[1], nofz)

        self.printInfo()
        

    ## simulate galaxies, returns numpy array of redshift, ra, dec
    def doSim(self):

        gal_data = np.zeros([self.ntotal, 3])
        for i in range(self.ntotal):

            gal_data[i,0],gal_data[i,1],gal_data[i,2] = simObject.drawObject(self)

        return gal_data


    ## Do simulation AND write to file
    # @param outfile    name of file to write to
    def writeSimToFile(self, outfile):

        f = open(outfile, 'w')

        gal_data = self.doSim()
        for i in range(self.ntotal):
            f.write(str(gal_data[i,0]) + '  ' + str(gal_data[i,1]) + '  ' + str(gal_data[i,2]) +'\n')

        f.close()


    ## Print basic info about simulation properties
    def printInfo(self):
        print 'Simulating background galaxies ...'
        simObject.printInfo(self)
        print ''

## Class to simulate clusters of galaxies
#
# Last modified 30 July 2014
class simulateClusters(simObject):


    ## The constructor
    # @param nclust        number of clusters to simulate
    # @param nofz
    # @param field_center  position in ra, dec (degrees) of field center
    # @param field_size    area of field in square degrees (assumed to be square)
    # @param zrange        simulate clusters within redshift range zrange=zlow,zhigh
    # @param richness      cluster richness (number of L* galaxies) range or constant value
    # @param mlim         magnitude limit of survey
    def __init__(self, nclust, nofz, field_center=(135.,30.), field_size=4., zrange=(1.,2.), \
                 richness=(10,201), mlim=26.5):
                 
        self.nclust = nclust
        self.mlim = mlim
        self.richness = richness
        self.isConstRich = False
        if (isinstance(richness,int)):
            self.isConstRich = True
    
        # init class that does the object position simulation
        ramin = field_center[0] - math.sqrt(field_size)/2.
        ramax = field_center[0] + math.sqrt(field_size)/2.
        decmin = field_center[1] - math.sqrt(field_size)/2.
        decmax = field_center[1] + math.sqrt(field_size)/2.
        simObject.__init__(self, ramin, ramax, decmin, decmax, zrange[0], zrange[1], nofz)
        
        self.printInfo()


    ## Simulate clusters, returns two numpy arrays (1) for clusters: cluster id, z, ra, dec, richness
    #  (2) for the cluster galaxies: cluster id, z, ra, dec
    def doSim(self):

        # array holding cluster data
        cl_data = np.zeros([self.nclust, 5])
        
        # loop over clusters
        for i in range(self.nclust):

            # draw cluster position and richness
            zcl,racl,deccl = simObject.drawObject(self)
            if (self.isConstRich):
                rcl = self.richness
            else:
                rcl = random.randint(self.richness[0],self.richness[1])
            
            # put in cluster data array
            cl_data[i,0],cl_data[i,1],cl_data[i,2],cl_data[i,3],cl_data[i,4] = i,zcl,racl,deccl,rcl
            
            # now draw cluster galaxies
            scg = simulateClusterGalaxies(rcl, (racl,deccl), zcl, self.mlim)
            g = scg.doSim(i) # returns array of cluster id, z, ra, dec
            
            # put into cluster galaxy data array
            if 'gal_data' in locals():
                gal_data = np.append(gal_data,g,axis=0)
            else:
                gal_data = g

        return cl_data,gal_data
        
    
    ## Do simulation AND write to file, one file contains list of clusters, other file contains
    #  the list of cluster member galaxies
    # @param outfileroot    root name of file to write to
    def writeSimToFile(self, outfileroot):

        # do simulation
        cl_data,gal_data = self.doSim()
        
        # write cluster data to file
        outfile = outfileroot + "_clusters.txt"
        f = open(outfile, 'w')
        for i in range(self.nclust):
            f.write(str(cl_data[i,0]) + '  ' + str(cl_data[i,1]) + '  ' + str(cl_data[i,2]) + '  ')
            f.write(str(cl_data[i,3]) + '  ' + str(cl_data[i,4]) + '\n')
        f.close()
        
        # write cluster member data to file
        outfile = outfileroot + "_clustermembers.txt"
        f = open(outfile, 'w')
        for i in range(len(gal_data)):
            f.write(str(gal_data[i,0]) + '  ' + str(gal_data[i,1]) + '  ' + str(gal_data[i,2]) + '  ')
            f.write(str(gal_data[i,3]) + '\n')
        f.close()
        
        
    ## Print basic info about simulation properties
    def printInfo(self):
        print 'Simulating galaxy clusters ...'
        if (self.isConstRich):
            print 'all with richness =',self.richness
        else:
            print 'with richness in the range',self.richness[0],'to',self.richness[1]
        simObject.printInfo(self)
        print ''


## Class to simulate cluster member galaxies
#
# Last modified 24 July 2014
class simulateClusterGalaxies(object):


    ## The constructor
    # @param richness     cluster richness (number of L* galaxies)
    # @param position     cluster position in ra and dec (degrees)
    # @param redshift     cluster redshift
    # @param mlim         magnitude limit of survey
    def __init__(self, richness, position, redshift, mlim=26.5):

        # cluster properties
        self.richness = richness
        self.ra_cluster = position[0]
        self.dec_cluster = position[1]
        self.z_cluster = redshift

        # angular diameter distance to cluster
        ccalc = cosmocalcs.cosmologyCalculator()
        ccalc.setEmissionRedShift(self.z_cluster)
        self.dA_cluster = ccalc.AngularDiameterDistanceMpc()
        
        # cluster density profile
        self.clusterProf = plummerProfile(redshift)
        self.rcut = self.clusterProf.rcut
        self.profmax = self.clusterProf.getProfile(0.)

        # get total number of L* cluster galaxies
        volCL = 4./3.*const.pi*self.rcut*self.rcut*self.rcut # volume of cluster optical richness definition
        clusterLF = galdist.schechterFunction()  # luminosity function
        Mstar = clusterLF.mstar          # mstar
        nLstar = self.richness*clusterLF.getphi(Mstar, redshift)*volCL # absolute number of L* galaxies
        self.nrcut =  int(round(self.richness*clusterLF.getnumberDensity(redshift, mlim)*volCL)) # in Mpc^-3 

        print 'Cluster richness =', self.richness ,'redshift =', self.z_cluster ,'position =', self.ra_cluster,
        print ',', self.dec_cluster,' a.d. distance =',self.dA_cluster,'Mpc'
        print 'Number of L* galaxies in cluster =', nLstar,' number of all in cluster =', self.nrcut,
        print ' out to', self.rcut 



    ## Draw cluster galaxy angular position (within rcut)
    def drawGalaxy(self):

        # draw galaxy radius
        rgal = -1.
        prof = -1. 
        while (True):
            rgal = random.uniform(0.,self.rcut)
            prof = random.uniform(0., self.profmax)

            p = self.clusterProf.getProfile(rgal)
            if (prof <= p):
                break
        #print 'Radius of cl gal',rgal
        # convert radius to angular position, randomized radial direction from cluster center
        phi = random.uniform(0.,2.*const.pi) # this is the angle in the triangle
        r_theta = rgal/self.dA_cluster       # this is: sqrt(dra^2 + ddec^2), ie r but in *radians*
        dra = r_theta*math.sin(phi)          # in radians
        ddec = r_theta*math.cos(phi)         # in radians
        ragal = self.ra_cluster - dra*const.rad2deg
        decgal = self.dec_cluster - ddec*const.rad2deg
       
        return ragal,decgal,r_theta*(180./const.pi)*60.


    ## simulate galaxies, returns numpy array of cluster id, z, ra, dec
    # @param cid    cluster id
    def doSim(self, cid=-1):

        minrad = 1000.
        maxrad = 0.
        meanrad = 0.
        gal_data = np.zeros([self.nrcut, 4])
        for i in range(self.nrcut):

            gal_data[i,0] = cid
            gal_data[i,1] = self.z_cluster
            gal_data[i,2],gal_data[i,3],tmp = self.drawGalaxy()
            meanrad+=tmp
            if (tmp>maxrad):
                maxrad = tmp
            if (tmp<minrad):
                minrad = tmp
        
        print 'Mean radius of cluster galaxy =',meanrad/self.nrcut,'max radius of cluster galaxy =',maxrad,
        print 'min radius of cluster galaxy =',minrad
        return gal_data


    ## Do simulation AND write to file
    # @param outfile    name of file to write to
    def writeSimToFile(self, outfile):

        f = open(outfile, 'w')

        gal_data = self.doSim()
        for i in range(self.nrcut):
            f.write(str(gal_data[i,0]) + '  ' + str(gal_data[i,1]) +'\n')

        f.close()
        
