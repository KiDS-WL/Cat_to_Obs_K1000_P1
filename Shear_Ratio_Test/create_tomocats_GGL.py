# ----------------------------------------------------------------
# File Name:           create_tomocats_GGL.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to select BOSS galaxies in fine tomographic
#                      bins for shear ratio test using data created by
#                      https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/GGL_LensCats/makelenscats.py
# ----------------------------------------------------------------
import sys
from astropy.io import fits 

# Read in user input to set the zmin zmax incatfile outcatfile
if len(sys.argv) <5: 
    print("Usage: %s zmin zmax catalogue.fits tomocat.fits" % sys.argv[0]) 
    sys.exit(1)
else:
    zmin = float(sys.argv[1]) 
    zmax = float(sys.argv[2]) 
    infile = sys.argv[3] 
    outfile = sys.argv[4]
    
#open the ldac catalogue using functions in ldac.py
#tests have shown ldac.py is much faster than using astropy

f = fits.open(infile)
#fits extension for GGL cats
iext=1

# read in the columns that we'll want 
# the spectroscopic redshift, weight, ra and dec

ALPHA_J2000=f[iext].data['ALPHA_J2000']
DELTA_J2000=f[iext].data['DELTA_J2000']
Z=f[iext].data['Z']
WEICOMP=f[iext].data['WEICOMP']

#This is the tomographic selection that we want to apply
ztomo=( (Z<zmax) & (Z>=zmin) )
ra_inbin=ALPHA_J2000[ztomo]
dec_inbin=DELTA_J2000[ztomo]
w_inbin=WEICOMP[ztomo]
z_inbin=Z[ztomo]

#Write out to output file - crucial that RA/DEC (in degrees) are double precision
#If you don't have that you round to a couple of arcsec for fields with ra > 100
hdulist = fits.BinTableHDU.from_columns(  
    [fits.Column(name='ALPHA_J2000', format='1D', unit='deg',array=ra_inbin),  
     fits.Column(name='DELTA_J2000', format='1D', unit='deg',array=dec_inbin),  
     fits.Column(name='WEICOMP', format='1E', array=w_inbin),
     fits.Column(name='Z', format='1E', array=z_inbin)])
hdulist.writeto(outfile, overwrite=True)
