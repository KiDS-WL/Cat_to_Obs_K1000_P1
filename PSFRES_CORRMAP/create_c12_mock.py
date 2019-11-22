import numpy as np
import sys
from astropy.io import fits
import ldac

md = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

##########################################################
# read in the c1 and c2 maps
hdu_list_c1 = fits.open(md+'/c1_map.fits')
c1_map = hdu_list_c1[0].data

hdu_list_c2 = fits.open(md+'/c2_map.fits')
c2_map = hdu_list_c2[0].data

# find the edges of the image from the exposure map
hdu_list_exp = fits.open(md+'/exposure_map.fits')
exp_map = hdu_list_exp[0].data

# and clip so it is zero for out of the image and 1 within
w_map = np.clip(exp_map,0,1)
w_map = w_map.astype(int)

# remove average c1/c2
dec1 = np.sum(w_map*c1_map)/np.sum(w_map)
dec2 = np.sum(w_map*c2_map)/np.sum(w_map)

c1_map = (c1_map - dec1)*w_map
c2_map = (c2_map - dec2)*w_map

###########################################################
#open the ldac catalogue using functions in ldac.py
ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the useful columns
KiDS_RA = ldac_table['ALPHA_J2000']
KiDS_Dec = ldac_table['DELTA_J2000']
Xpos_in = ldac_table['Xpos_THELI']
Ypos_in = ldac_table['Ypos_THELI']
weight_in = ldac_table['recal_weight']

Nobj = np.shape(Xpos_in)[0]

XYpos = np.vstack((Xpos_in,Ypos_in))
XYpos_short = XYpos/10.
XYpos_int = XYpos_short.astype(int)

# mock e1 and e1 columns
c1 = np.zeros(Nobj)
c2 = np.zeros(Nobj)

for i in range(Nobj):
        c1[i] = c1_map[XYpos_int[0,i],XYpos_int[1,i]]
        c2[i] = c2_map[XYpos_int[0,i],XYpos_int[1,i]]

ldac_table['c1'] = c1
ldac_table['c2'] = c2

ldac_table.saveas(outfile)
