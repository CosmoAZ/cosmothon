import numpy as np


def group(pgroup, plist, b):
    """ 
        pgroup: particles in group list of (id,x_i,x_j,x_k, ... x_n)
        plist: list of particles not in group list of (id,x_i,x_j,x_k, ... x_n)
        
        where x_i to x_n are the coordinates of the particle in n dimensional space
        
        b: linking length
    """
    
    # record initial lengths of lists
    ngroup = len(pgroup)
    nlist = len(plist)
    
    # id of everything in the group
    idg = np.array([i[0] for i in pgroup])  # O(ng)
    
    # x,y positions of everything in the group
    #xg = np.array([i[1] for i in pgroup])   # O(ng)
    #yg = np.array([i[2] for i in pgroup])   # O(ng)
    pg = np.array([i[1:] for i in pgroup])

    # loop over particles NOT in the group
    for p in plist:
    
        # particle properties
        idp = p[0] # +1
        #xp = p[1]  # +1
        #yp = p[2]  # +1
        pp = np.array(p[1:])

        # distance between this particle and all in the group already
        dr = np.sqrt(np.sum((pg-pp)*(pg-pp),axis=1))
        #np.sqrt((xg-xp)*(xg-xp) + (yg-yp)*(yg-yp)) # +7
        
        # indices of the particles with distance < b from a particle in the group
        id_in_group = np.where(dr<b)[0] # +f(nl)
        
        # if there is at least one particle in the group within b of this particle 
        if (len(id_in_group)>0):
            
            # add this particle to the group
            pgroup.append(p)  # +1
            # remove matching particle from list
            plist.remove(p)   # +1
                

    # check if nothing changed
    if (abs(len(pgroup)-ngroup)<0.1):
        return pgroup, plist
    else:
        return group(pgroup, plist, b)
    
    
