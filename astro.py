# useful astronomy methods

import math
#import const

def convert_to_deg(rah,ram,ras,decd,decm,decs):
    """ RA and Dec. coordinate conversion from RA in hours,mins,secs and Dec in
        deg,arcmins,arcsecs to decimal degrees
		
		(num,num,num,num,num,num) -> (num,num)
		
    """
    
    ra_deg = 15.*(rah + ram/60. + ras/3600.)
    if decd>=0.:
        dec_deg = decd + decm/60. + decs/3600.
    else:
        dec_deg = decd - decm/60. - decs/3600.
    
    return ra_deg,dec_deg
