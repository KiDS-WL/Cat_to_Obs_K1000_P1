import astropy.io.fits as fits
import numpy as np
import sys
import string
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.wcs as pywcs
import astropy

#This script takes the 2dFLenS mask (in low res ascii format) and pastes it
#onto the high res KiDS mask

#input data:
#we have the high rest KiDS mask (fits) and the 2dFLenS mask (currently ascii)
#Download KiDS mask from http://cuillin.roe.ac.uk/~cech/KiDS/KiDS-1000/MASKS_and_effective_area/
#It's to big to add to the repo directly

KiDS_mask_fits='KiDS-S_K1000_footprint_16bit.fits'
twodFLenS_mask_asc='mask_kidss.dat'

# read in the 2dFLenS ascii mask
maskdata = np.loadtxt(twodFLenS_mask_asc)
ra=maskdata[:,0]
dec=maskdata[:,1]
maskflag=maskdata[:,2]

#the 2dFLenS mask ascii file has ra spaced by 7 arcsec (0.11666666666666667 degrees)
#the 2dFLenS mask ascii file has dec spaced by 0.1 degrees
#ra starts at 330.058
#dec starts at -35.950

#unfortunately it ends at ra=52.442 before the KiDS-S footprint ends which adds some complexity
#when we put the two masks together 
#easiest route is to add in 22 extra RA points to extend the mask to RA=55
#Dec runs from -26.050 to -35.950 = 100 bins

ra_add=np.zeros(22*100)
dec_add=np.zeros(22*100)
maskflag_add=np.zeros(22*100)  # no 2dFLenS imaging here

k=0
for i in range(22):
    for j in range(100):
        ra_add[k]= 52.442 +(i+1)*7.0/60.0
        dec_add[k] = -35.950 + j*0.1
        k = k + 1
        
ra = np.append(ra,ra_add)
dec = np.append(dec,dec_add)
maskflag = np.append(maskflag, maskflag_add)

#place the 2dFLenS mask into an array
#include this trick to deal with zero crossing for KiDS-S
ra[ra < 180.] += 360.
# form an integer grid
imask_x = np.rint((ra-330.058)*60.0/7.0)
imask_y = np.rint((dec+35.950)/0.1)
nmaskpix_x = np.max(imask_x).astype(int)+1
nmaskpix_y = np.max(imask_y).astype(int)+1

#now we know how our grid corresponds to ra/dec lets fill up the array
twodFLenS_mask=np.ones((nmaskpix_x,nmaskpix_y),dtype=int)
for i in range(len(imask_x)):
    if maskflag[i]==1:  # imaged by 2dFLenS
        twodFLenS_mask[np.int(imask_x[i]),np.int(imask_y[i])]=0
    else:  #not imaged by 2dFLenS
        twodFLenS_mask[np.int(imask_x[i]),np.int(imask_y[i])]=1
        
#Time to append this information to the KiDS mask which is on a finer resolution

#First read in the KiDS mask and the wcs
KiDSmaskimage = fits.open(KiDS_mask_fits) # axis flipped!
KiDSmaskdata = KiDSmaskimage[0].data
wcs = WCS(KiDS_mask_fits)

# calculate the ra and dec of every x,y position in the KiDS fits which are not outside the wcs cut
# dimensions of the KiDS fits mask
naxis_x=72000
naxis_y=12000

#set up the new joint mask - starting with the KiDS mask
#also set up an empty mask of 1's with the same dimensions to make the 2dFLenS-only high-res version

jointmask = np.copy(KiDSmaskdata)
twodFLenS_mask_hires = np.ones_like(KiDSmaskdata) 

#you may only be interested in the pixels that are in the KiDS footprint
#infoot=np.where(KiDSmaskdata<16384)
#or you may want to see the full extent of 2dFLenS which extends beyond KiDS, in which case
infoot=np.where(KiDSmaskdata>=0)

# given the wcs, get the ra/dec for these pixels
coords=astropy.wcs.utils.pixel_to_skycoord(infoot[1], infoot[0], wcs, origin=0, mode='all', cls=None)

# now find which pixel in the 2dFLenS mask this ra/dec corresponds to
# need to include the same zero crossing trick
rakids=np.array(coords.ra)
rakids[rakids < 180.] += 360.
    
imask_x = np.rint((rakids-330.058)*60.0/7.0)
imask_y = np.rint((np.array(coords.dec)+35.950)/0.1)

# the KiDS mask extends beyond the 2dFLenS mask but we need to only look at pixels where
# the twodFLenS_mask is defined

gdtwodflenspix=((imask_x >= 0) & (imask_x<nmaskpix_x) & (imask_y >= 0) & (imask_y<nmaskpix_y))

# We've run out of bitmask values for this joint mask
# add in the 2dFLenS mask using a value "8" which is also used for the manual masking of globular clusters over a very small area
jointmask[infoot[0][gdtwodflenspix], infoot[1][gdtwodflenspix]] = KiDSmaskdata[infoot[0][gdtwodflenspix], infoot[1][gdtwodflenspix]] + twodFLenS_mask[imask_x.astype(int)[gdtwodflenspix],imask_y.astype(int)[gdtwodflenspix]]*8

#By adding the 2dFLenS 8-bit mask there will be a handful of previously manually masked 8-bit
#pixels that have maxed out the bitpix flagging in every mask 
# i.e not in 2dFLenS, and manually masked and everything else gone wrong aswell!
#we need to remove those now negative points - it's around ~100 pixels

jointmask[jointmask < 0] = 32767

# lets also make our high-res 2dFLenS mask
twodFLenS_mask_hires[infoot[0][gdtwodflenspix], infoot[1][gdtwodflenspix]] = twodFLenS_mask[imask_x.astype(int)[gdtwodflenspix],imask_y.astype(int)[gdtwodflenspix]]

#finally write out the joint mask
KiDSmaskimage[0].data = jointmask
KiDSmaskimage.writeto('KiDS-S_K1000_w_2dFLenS_footprint_full.fits',overwrite=True)

#and write out the 2dFLenS high res mask
KiDSmaskimage[0].data = twodFLenS_mask_hires
KiDSmaskimage.writeto('KiDS-S_2dFLenS_only_footprint.fits',overwrite=True)

#This output is now stored in http://cuillin.roe.ac.uk/~cech/KiDS/KiDS-1000/MASKS_and_effective_area/
#The effective unmasked area is calculated when this is converted to healpix in healpix_mosaic_mask.py
