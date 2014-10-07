""" Simple python implementation of "friends-of-friends" algorithm to find clustered objects
    (Huchra & Geller 1982)

"""

import numpy as np


def group(pgroup, plist, b):
    """ 
        pgroup: particles already in group, list of tuples in form (id,x_i,x_j,x_k, ... x_n)
        plist:  particles not in yet group, list of tuples in form (id,x_i,x_j,x_k, ... x_n)

        where
        id is the particle id
        x_i to x_n are the coordinates of the particle in n dimensional space
        
        b: linking length (how close a particle needs to be to be a "friend")
    """
    
    # record initial lengths of lists
    ngroup = len(pgroup)
    nlist = len(plist)
    
    # id of everything in the group
    idg = np.array([i[0] for i in pgroup]) 
    
    # positions of everything in the group
    pg = np.array([i[1:] for i in pgroup])

    # loop over particles NOT in the group
    for p in plist:
    
        # particle properties
        idp = p[0] # particle id
        pp = np.array(p[1:]) # particle position

        # distance between this particle and all in the group already
        dr = np.sqrt(np.sum((pg-pp)*(pg-pp),axis=1))

        # indices of the particles with distance < b from a particle in the group
        id_in_group = np.where(dr<b)[0] 
        
        # if there is at least one particle in the group within b of this particle 
        if (len(id_in_group)>0):
            
            # add this particle to the group
            pgroup.append(p)  
            # remove matching particle from list
            plist.remove(p)   
                

    # check if nothing changed
    if (abs(len(pgroup)-ngroup)<0.1):
        return pgroup, plist
    else:
        return group(pgroup, plist, b)
    
    
