""" Demo of friends-of-friends algorithm

    outputs to file:
    particles.txt: x,y positions of all particles to group
    groups.txt: each row contains list of x y of all particles in a group, one row per group

"""
import sys, getopt
import math
import numpy as np
import random
import time

import fof


    
def main(argv):
    
    ## create list of particles with x,y positions
    ## they will be in two clusters with x,y centers below
    
    cluster1 = (4.,5.)
    cluster2 = (1.,2.)

    # number of particles in both clusters
    npoints = 500 
    # list of particles
    plist = []
    # clusters have width 0.5, Gaussian distributed

    # cluster 1 
    for i in xrange(0,npoints/2):

        x = random.gauss(cluster1[0],0.5)
        y = random.gauss(cluster1[1],0.5)
        idg = i
        plist.append((idg,x,y))

    # cluster 2
    for i in xrange(npoints/2,npoints):
        x = random.gauss(cluster2[0],0.5)
        y = random.gauss(cluster2[1],0.5)
        idg = i
        plist.append((idg,x,y))

        
    ## add some random background

    # number of background particles
    nback = 1000
    xmin = min([i[1] for i in plist])
    xmax = max([i[1] for i in plist])
    ymin = min([i[2] for i in plist])
    ymax = max([i[2] for i in plist])
    print xmin,'< x <',xmax,ymin,'< y <',ymax
    for i in xrange(npoints,npoints+nback):
        x = random.uniform(xmin,xmax)
        y = random.uniform(ymin,ymax)
        idg = i
        plist.append((idg,x,y))

        
    ## write full distribution to a file
    f = open('particles.txt','w')
    for p in plist:
        f.write(str(p[0]) + '  ' + str(p[1]) + '  ' + str(p[2]) + '\n')
    f.close()


    ## Run friends-of-friends
       
    # linking length
    b = 0.125
    groups = []
    
    i = 1;
    t1 = time.time()
    ng = 0
    while len(plist)>0:
        print 'group',i,
        
        # initialize  
        pgroup = [plist[0]]
        plist.remove(pgroup[0])
        
        g, plist = fof.group(pgroup, plist, b)
        ng += len(g)
        print len(g),'galaxies in the group',ng,'galaxies grouped so far'
        groups.append(g)
        i+=1

    print 'Found',len(groups),'groups from',npoints+nback,'galaxies'
    t2 = time.time()
    print 'Time taken',t2-t1,'s'
    
    
    ## write groups to a file
    outfile = 'grouped.txt'
    f = open(outfile,'w')
    ii = 1
    for gr in groups:
    
        print len(gr),'in group',ii
        for p in gr:
            f.write(str(p[1]) + '  ' + str(p[2]) + '  ')
        f.write('\n')
        ii+=1
    f.close()


# Run program
if __name__ == '__main__':

    main(sys.argv[1:])
    
    
