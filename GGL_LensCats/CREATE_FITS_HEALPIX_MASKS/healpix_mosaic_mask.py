from astropy.io import fits
from astropy.wcs import WCS
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from astropy import units as u
from astropy.coordinates import SkyCoord

#========================================
#2dFLenS_ascmask_to_fits.py creates fits masks at high res
#in this script we convert them into healpix masks
#and calculate the effective area of the overlap
#========================================

#1: Read in the fitsmask and set the value you wish to mask out

#Either read 2dFLenS flat mask
#twodFfile='/home/cech/KiDSLenS/THELI_catalogues/MOSAIC_MASK/HEALPIX_MASK/KiDS-S_2dFLenS_only_footprint.fits'

# or the joint 2dFLenS-KiDS mask
twodFfile='/home/cech/KiDSLenS/THELI_catalogues/MOSAIC_MASK/HEALPIX_MASK/KiDS-S_K1000_w_2dFLenS_footprint_full.fits'

#Set the bitmask that we are applying
#info here:http://lensingkids.strw.leidenuniv.nl/doku.php?id=kids-1000#mask_values_and_meanings

#For the joint KiDS-2dFLenS mask this will give the overlap
#of sources that have 9-band photo-zs and 2dFLenS coverage
#bitmask=0x6FFC

#You might however just want to know where the gri information is
#for the 2dFLenS re-weighting to BOSS - for this you want
#bitmask=0x681C

#Maybe you just want to know if VST has ever pointed at this
#patch of sky = i.e the KiDS extend (wcs mask - 16384 + 8 for the 2dFLenS coverage)
bitmask=0x4008

#Or maybe you are looking at the 2dFLenS only footprint
#bitmask=0x1

#Read in the mask file
hdulist = fits.open(twodFfile)
wcs_data = WCS(hdulist[0].header)
image_data = hdulist[0].data
ifilter=np.logical_not(np.array(image_data & bitmask, dtype=bool))

#======================================
#2- Link pixels with RA,Dec

x,y=np.arange(1,hdulist[0].header['NAXIS1']+1),np.arange(1,hdulist[0].header['NAXIS2']+1)
xx, yy = np.meshgrid(x, y)
ra,dec=wcs_data.wcs_pix2world(xx[ifilter],yy[ifilter],1)

#Use this byproduct of the healpix production to ask what is the
#effective area of the unmasked data?
#Count the pixels that are unmasked and multiply by the pixel scale
#which is 0.1x0.1 arcmin^2 (see CD1_1 and CD1_2 in header of fitsmask)

unmasked_eff_area=len(ra)*(0.1*0.1)/(60.0*60.0)

#bitmask=1
#print ('2dFLenS-SGP effective area = %f sq. degrees'%(unmasked_eff_area))
# 9-band bitmask
#print ('2dFLenS-KiDS-S 9-band effective area = %f sq. degrees'%(unmasked_eff_area))
# gri bitmask
#print ('2dFLenS-KiDS-S gri effective area = %f sq. degrees'%(unmasked_eff_area))
# wcs bitmask
print ('2dFLenS-KiDS-S wcs effective area = %f sq. degrees'%(unmasked_eff_area))

radec=np.stack((ra,dec)).T

theta = np.deg2rad(90.0 - radec[:, 1])
phi = np.deg2rad(radec[:, 0])

#===========================================

#3 - Create a "nside" healpix template map and write to file
nside=4096
nest=False
mask_hppix = hp.ang2pix(nside, theta=theta, phi=phi, nest=nest)
npix = hp.nside2npix(nside)
weights=None
countmap = np.bincount(mask_hppix, weights=weights, minlength=npix)

# Write out the 2dFLenS only footprint
#hp.write_map('2dFLenS_SGP_healpix.fits', countmap,coord='G',overwrite=True)
#hp.mollview(countmap,coord='C',xsize=20000, title='2dFLenS-SGP')
#plt.savefig('2dFLenS_SGP_healpix.png')

#Or write out the joint lensing-2dFLenS 9-band overlap
#hp.write_map('K1000_w_2dFLenS_9_band_overlap_healpix.fits', countmap,coord='G',overwrite=True)
#hp.mollview(countmap,coord='C',xsize=20000, title='KiDS-2dFLenS 9-band overlap')
#plt.savefig('K1000_w_2dFLenS_9_band_overlap.png')

#Or write out the joint lensing-2dFLenS gri overlap
#hp.write_map('K1000_w_2dFLenS_gri_overlap_healpix.fits', countmap,coord='G',overwrite=True)
#hp.mollview(countmap,coord='C',xsize=20000, title='KiDS-2dFLenS gri overlap')
#plt.savefig('K1000_w_2dFLenS_gri_overlap.png')

#Or write out the joint lensing-2dFLenS KiDS extemt overlap
hp.write_map('K1000_w_2dFLenS_wcs_overlap_healpix.fits', countmap,coord='G',overwrite=True)
hp.mollview(countmap,coord='C',xsize=1000, title='KiDS-2dFLenS wcs overlap')
plt.savefig('K1000_w_2dFLenS_wcs_overlap.png')
