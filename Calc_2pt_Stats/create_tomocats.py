# ----------------------------------------------------------------
# File Name:           create_tomocats.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to select galaxies in a tomographic
#                      bin and write out a new file with only the columns we need 
#                      helps with I/O to only have small cats
#                      if CCORR=true then we apply a c-correction
#                      script will need to change if keywords in KIDS cats are updated
#                      For metacal we'll also need to pass through m
# ----------------------------------------------------------------
import sys
import ldac
import numpy as np
from astropy.io import fits 

# Read in user input to set the zmin zmax incatfile outcatfile CCORR?
if len(sys.argv) <7: 
    print "Usage: %s zmin zmax catalogue.fits tomocat.fits blind_ID ccorr(true or false)" % sys.argv[0] 
    sys.exit(1)
else:
    zmin = float(sys.argv[1]) 
    zmax = float(sys.argv[2]) 
    infile = sys.argv[3] 
    outfile = sys.argv[4]
    blind = sys.argv[5]
    ccorr = sys.argv[6]
    
#open the ldac catalogue using functions in ldac.py
#tests have shown ldac.py is much faster than using astropy
ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the blinded ellipticity columns that we'll want 
# and the photometric redshift, ra and dec
e1colname='e1_'+blind
e2colname='e2_'+blind
wtcolname='weight_'+blind    
e1=ldac_table[e1colname]
e2=ldac_table[e2colname]
weight=ldac_table[wtcolname]
Z_B=ldac_table['Z_B']

ALPHA_J2000=ldac_table['ALPHA_J2000']
DELTA_J2000=ldac_table['DELTA_J2000']

#This is the tomographic selection that we want to apply
ztomo=( (Z_B<zmax) & (Z_B>=zmin))
e1_inbin=e1[ztomo]
e2_inbin=e2[ztomo]
w_inbin=weight[ztomo]
ra_inbin=ALPHA_J2000[ztomo]
dec_inbin=DELTA_J2000[ztomo]

# THIS WILL NEED TO BE UBDATED FOR METACAL TO ALSO INCLUDE THE M_CORRECTION
if (ccorr=='true'):  
    # weighted mean   
    c1=np.average(e1_inbin,weights=w_inbin)
    c2=np.average(e2_inbin,weights=w_inbin)
    # weighted sd
    sumosq1=np.average(e1_inbin*e1_inbin,weights=w_inbin)
    sumosq2=np.average(e2_inbin*e2_inbin,weights=w_inbin)
    errc1= np.sqrt(sumosq1 - c1*c1)
    errc2= np.sqrt(sumosq2 - c2*c2)

    print zmin, zmax, c1, errc1, c2, errc2  
else:
    c1 = 0
    c2 = 0

#Apply correction
e1_corr = e1_inbin - c1
e2_corr = e2_inbin - c2

#Write out to output file - crucial that RA/DEC (in degrees) are double precision
#If you don't have that you round to a couple of arcsec for fields with ra > 100
hdulist = fits.BinTableHDU.from_columns(  
    [fits.Column(name='ALPHA_J2000', format='1D', unit='deg',array=ra_inbin),  
     fits.Column(name='DELTA_J2000', format='1D', unit='deg',array=dec_inbin),  
     fits.Column(name='e1', format='1E', array=e1_corr),  
     fits.Column(name='e2', format='1E', array=e2_corr),  
     fits.Column(name='weight', format='1E', array=w_inbin)])
hdulist.writeto(outfile, overwrite=True)
