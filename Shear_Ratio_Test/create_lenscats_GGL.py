# ----------------------------------------------------------------
# File Name:           create_lenscats_GGL.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to select BOSS galaxies in fine tomographic
#                      bins for shear ratio test using data created by
#                      https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/GGL_LensCats/makelenscats.py
# ----------------------------------------------------------------
import sys
import numpy as np
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

# Whether input is FITS or ascii, the FITS extension and keywords for
# ra,dec,weight,redshift, depend on input catalogue. Que if statements:

if "MICE2" in infile:
    Read_FITS=True            # Input is a FITS file
    iext=1
    ra_keyword='ra_gal'
    dec_keyword='dec_gal'
    w_keyword='recal_weight'
    z_keyword='z_cgal_v'
    if "_on_" in infile:       # Magnification is on.
        ra_keyword += '_mag'
        dec_keyword += '_mag'

elif "BOSS" in infile:
    Read_FITS=True              # Input is a FITS file
    iext=1                      #fits extension

    ra_keyword='ALPHA_J2000'
    dec_keyword='DELTA_J2000'
    w_keyword='WEICOMP'
    z_keyword='Z'

elif "GAMA" in infile:
    Read_FITS=False             # Input is an ascii file
                                # Assume columns ordered (RA,Dec,z,w)
        
else:
    print("This code only accepts infiles containing 'BOSS', 'GAMA', 'MICE2' keywords.")
    print("This infile ",infile," contains none of these, so I'm exiting.")
    sys.exit(1)

if Read_FITS:
    # Read in the columns that we'll want 
    # the spectroscopic redshift, weight, ra and dec
    f = fits.open(infile)
    ALPHA_J2000=f[iext].data[ra_keyword]
    DELTA_J2000=f[iext].data[dec_keyword]
    Z=f[iext].data[z_keyword]
    WEICOMP=f[iext].data[w_keyword]
else:
    # Read in an ascii file
    ALPHA_J2000,DELTA_J2000,Z,WEICOMP=np.loadtxt(infile, usecols=(0,1,2,3), unpack=True)

    
# HACK. UNCOMMENT THIS TO ONLY PICK THE MASSIVE LENSES.
#heavies = np.where( f[iext].data['lmstellar'] > 11.5 )
#ALPHA_J2000 = ALPHA_J2000[heavies]
#DELTA_J2000=DELTA_J2000[heavies]
#Z=Z[heavies]
#WEICOMP=np.ones_like(Z) #WEICOMP[heavies]
#print(WEICOMP)
# END OF HACK.

if "MICE2" in infile and "random" in outfile:
    print("We are making MICE randoms. The outfile is %s"%outfile)
    f = 1  # Scale ngals*f random points
    # Shuffle the positions of the gals in the MICE mocks...
    # ...replace the RA,Dec's with random values
    np.random.seed(42)
    ra_rand = np.random.uniform( ALPHA_J2000.min(), ALPHA_J2000.max(), f*len(ALPHA_J2000) )
    # To avoid Dec clustering at the poles, uniformly generate in sin(Dec)
    dec_rad = DELTA_J2000 * np.pi / 180. # Dec in radians
    np.random.seed(43)
    dec_rand = np.arcsin( np.random.uniform( np.sin(dec_rad.min()), np.sin(dec_rad.max()), f*len(dec_rad) ) ) 
    DELTA_J2000 = dec_rand * 180. / np.pi # convert back to deg
    #DELTA_J2000 = np.random.uniform( DELTA_J2000.min(), DELTA_J2000.max(), f*len(DELTA_J2000) )
    ALPHA_J2000 = np.copy( ra_rand )
    WEICOMP = np.ones( f*len(WEICOMP) )


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
