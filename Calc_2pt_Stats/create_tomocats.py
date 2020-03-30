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
    print("Usage: %s zmin zmax catalogue.fits tomocat.fits blind_ID ccorr(true or false)" % sys.argv[0])
    sys.exit(1)
else:
    zmin = float(sys.argv[1])
    zmax = float(sys.argv[2])
    infile = sys.argv[3]
    outfile = sys.argv[4]
    blind = sys.argv[5]
    ccorr = sys.argv[6]
    
if len(sys.argv) >7:     # optional commands to carry through the SOM Flag of choice
    flag_SOM = sys.argv[7]
else:
    flag_SOM = False

# Define the bootstrap error function to calculate the error on the c-terms
def Bootstrap_Error(nboot, samples, weights):
	N = len(samples)
	bt_samples = np.zeros(nboot)		 		# Will store mean of nboot resamples
	for i in range(nboot):
		idx = np.random.randint(0,N,N)			# Picks N random indicies with replacement
		bt_samples[i] = np.sum( weights[idx]*samples[idx] ) / np.sum( weights[idx])
	return np.std(bt_samples)

#open the ldac catalogue using functions in ldac.py
#tests have shown ldac.py is much faster than using astropy
ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the blinded ellipticity columns that we'll want
# and the photometric redshift, ra and dec
# using phase 1 autocal column names
e1colname='autocal_e1_'+blind
e2colname='autocal_e2_'+blind
wtcolname='recal_weight_'+blind 

#phase 0 column names  
#e1colname='e1_'+blind
#e2colname='e2_'+blind
#wtcolname='weight_'+blind

e1=ldac_table[e1colname]
e2=ldac_table[e2colname]
weight=ldac_table[wtcolname]
Z_B=ldac_table['Z_B']

# Lets also pass through the PSF ellipticity for star-gal-xcorr 
PSF_e1=ldac_table['PSF_e1']
PSF_e2=ldac_table['PSF_e2']

ALPHA_J2000=ldac_table['ALPHA_J2000']
DELTA_J2000=ldac_table['DELTA_J2000']

# If a SOM Flag has been selected, them we need to apply the SOM Flag in the
# tomographic selection
if flag_SOM:
    flag_SOM_name=flag_SOM+'_'+blind
    FLAG_SOM=ldac_table[flag_SOM_name]
    ztomo=( (Z_B<=zmax) & (Z_B>zmin) & (FLAG_SOM > 0))
else:
    ztomo=( (Z_B<=zmax) & (Z_B>zmin))
    
#Apply the tomographic/SOM selection

e1_inbin=e1[ztomo]

e2_inbin=e2[ztomo]
ra_inbin=ALPHA_J2000[ztomo]
dec_inbin=DELTA_J2000[ztomo]
PSF_e1_inbin=PSF_e1[ztomo]
PSF_e2_inbin=PSF_e2[ztomo]
w_inbin=weight[ztomo]

#carry through the square of the weight for
#Npair calculation hack with Treecorr
wsq_inbin=weight[ztomo]*weight[ztomo]

nboot = 300

# THIS WOULD NEED TO BE UPDATED FOR METACAL TO ALSO INCLUDE THE M_CORRECTION
if (ccorr=='true'):  
    # weighted mean   
    c1=np.average(e1_inbin,weights=w_inbin)
    c2=np.average(e2_inbin,weights=w_inbin)

    # Bootstrap error on the mean
    errc1=Bootstrap_Error(nboot, e1_inbin, w_inbin)
    errc2=Bootstrap_Error(nboot, e2_inbin, w_inbin)

    # weighted sd
    #sumosq1=np.average(e1_inbin*e1_inbin,weights=w_inbin)
    #sumosq2=np.average(e2_inbin*e2_inbin,weights=w_inbin)
    #errc1= np.sqrt(sumosq1 - c1*c1)
    #errc2= np.sqrt(sumosq2 - c2*c2)

    print(zmin, zmax, c1, errc1, c2, errc2)
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
     fits.Column(name='PSF_e1', format='1E', array=PSF_e1_inbin),
     fits.Column(name='PSF_e2', format='1E', array=PSF_e2_inbin),
     fits.Column(name='weight', format='1E', array=w_inbin),
     fits.Column(name='weightsq', format='1E', array=wsq_inbin)])
hdulist.writeto(outfile, overwrite=True)

#ascii output for athena
#np.savetxt('test.asc',np.transpose([ra_inbin,dec_inbin,e1_corr,e2_corr,w_inbin]),header='ALPHA_J2000 DELTA_J2000 e1 e2 weight',fmt='%.12e')  
