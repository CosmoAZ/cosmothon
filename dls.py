"""
    Functions specific to the Deep Lens Survey


"""
import wcslib
#import MySQLdb
import pyfits
import astro
import math


# pixel scale of mosaic camera
pixsc = 0.257 

# wcs info
fileF1 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F1wcs.fits')
hdrF1 = fileF1[0].header
fileF2 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F2wcs.fits')
hdrF2 = fileF2[0].header
fileF3 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F3wcs.fits')
hdrF3 = fileF3[0].header
fileF4 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F4wcs.fits')
hdrF4 = fileF4[0].header
fileF5 = pyfits.open('/home/alex/DeepLensSurvey/wcs/F5wcs.fits')
hdrF5 = fileF5[0].header

# Projection parameters from header
wcsF1 = wcslib.WcsProjection(hdrF1)
wcsF2 = wcslib.WcsProjection(hdrF2)
wcsF3 = wcslib.WcsProjection(hdrF3)
wcsF4 = wcslib.WcsProjection(hdrF4)
wcsF5 = wcslib.WcsProjection(hdrF5)


def returnxy(field_number,ra,dec):
    """ for a given DLS field number, object ra, dec and the WCS transformation
        systems for each DLS field, returns the pixel coordinates

        field_number is in the format 'F#'

        ra,dec are in degrees

    """

    field_list = ['F1','F2','F3','F4','F5']
    if field_number in field_list:

        print 'Converting ra dec via WCS transformation for field',field_number
        if field_number == 'F1':
            x,y = wcsF1.toPixel(ra, dec, rnd=False)
        elif field_number == 'F2':
            x,y = wcsF2.toPixel(ra, dec, rnd=False)
        elif field_number == 'F3':
            x,y = wcsF3.toPixel(ra, dec, rnd=False)
        elif field_number == 'F4':
            x,y = wcsF4.toPixel(ra, dec, rnd=False)
        else:
            x,y = wcsF5.toPixel(ra, dec, rnd=False)

    else:
        print "Error: field",field_number,"not understood"
            
    return x,y


def returnFieldLimits(field):
    """ Return the min and max ra,dec in degrees for the given DLS field
    """

    field_list = ['F1','F2','F3','F4','F5']
    if field in field_list:

        if field == 'F1':
            rah,ram,ras=0.,53.,25.3 	
            decd,decm,decs=12.,33.,55. 
            ra,dec = astro.convert_to_deg(rah, ram, ras, decd, decm, decs)
        elif field == 'F2':
            rah,ram,ras=9.,19.,32.4
            decd,decm,decs=30.,0.,0.
            ra,dec = astro.convert_to_deg(rah, ram, ras, decd, decm, decs)
        elif field == 'F3':
            rah,ram,ras=5.,20.,0.
            decd,decm,decs=-49.,0.,0. 
            ra,dec = astro.convert_to_deg(rah, ram, ras, decd, decm, decs)
        elif field == 'F4':
            rah,ram,ras=10.,52.,0.
            decd,decm,decs=-5.,0.,0. 
            ra,dec = astro.convert_to_deg(rah, ram, ras, decd, decm, decs)
        else:
            rah,ram,ras=13.,59.,20.
            decd,decm,decs=-11.,3.,0.
            ra,dec = astro.convert_to_deg(rah, ram, ras, decd, decm, decs)

    else:
        print "Error: field",field,"not understood"

    ramin = ra-1.
    ramax = ra+1.
    decmin = dec-1.
    decmax = dec+1.
    lims = [ramin,ramax,decmin,decmax]
    return lims
    
def checkInDLS(ra,dec):
    """ for a given ra,dec checks whether this position lies in any DLS field 
    """

    isInDLS = False
    field_list = ['F1','F2','F3','F4','F5']
    for field in field_list:
        #print 'Checking field ',field
        [ramin,ramax,decmin,decmax] = returnFieldLimits(field)
        if (ra>=ramin and ra<=ramax and dec>=decmin and dec<=decmax):
            isInDLS = True
            break
            
    return isInDLS

def findField(ra,dec,cursor):
    """ for a given ra, dec return the DLS field they are located in
    
        field number returned is in the format 'F#'

        ra,dec are in degrees

    """

    field = 'FF'
    if (not checkInDLS(ra,dec)):
        return field

    getSubfield = 'select DLSref.fNearestSubfield(' + str(ra) + ',' + str(dec) + ')' 
    cursor.execute(getSubfield)
    results = cursor.fetchall()
    field = results[0][0][:2]
    return field
    

def convertDecimalToBinary(x):
    """ convert decimal number x into binary
    """
    
    if (x<1):
        binary_value = [0]
        return binary_value
    
    binary_holder = []
    while (math.floor(x/2.) > 0):
        binary_holder.append(x%2)
        x /= 2.
    
    if (len(binary_holder)<1):
        binary_value = [1]
        return binary_value
    else:
        binary_holder.append(1.)
    
    binary_value = []
    for binary in reversed(binary_holder):
        binary_value.append(int(binary))
    
    return binary_value
