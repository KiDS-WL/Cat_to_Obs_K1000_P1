#
import sys
import numpy as np
import ldac

"""This script takes an input LDAC catalogue
"""
#

# Read in user input
if len(sys.argv) <1: 
    print "Usage: %s InputCat field " % sys.argv[0] 
    sys.exit(1)
else:
    field=sys.argv[1]
    infile="/disk09/KIDS/KIDSCOLLAB_V1.0.0/%s/r_SDSS/colourcat_V1.0.0A/%s_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_sel.cat" %(field,field)


#open the ldac catalogue using functions in ldac.py
ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the m and ellipticity columns to be binned
e1colname='raw_e1'
e2colname='raw_e2'
m1colname='metacal_m1'
m2colname='metacal_m2'
magname='catmag'

e1=ldac_table[e1colname]
e2=ldac_table[e2colname]
m1=ldac_table[m1colname]
m2=ldac_table[m2colname]
mag=ldac_table[magname]

ngals_all = len(e1)
#print 'Read in ',ngals_all,' raw ellipticities and m1/m2 from ',infile

nbins=20
magmin=20
magmax=25

ngal_bin,magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),density=False)
m1_bin,magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),weights=m1,density=False)
m2_bin,magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),weights=m2,density=False)
e1_sq_bin, magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),weights=e1*e1,density=False)
e2_sq_bin, magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),weights=e2*e2,density=False)
avemag_bin,magbin = np.histogram(mag,bins=nbins,range=(magmin,magmax),weights=mag,density=False)


m1_bin = m1_bin/ngal_bin
m2_bin = m2_bin/ngal_bin
sigma_e1_bin = np.sqrt(e1_sq_bin/ngal_bin)
sigma_e2_bin = np.sqrt(e2_sq_bin/ngal_bin)
avemag_bin = avemag_bin/ngal_bin

print '# bin_int avemag_bin m1_bin m2_bin sigma_e1_bin sigma_e2_bin ngal_bin'
for i in range(nbins):
    print i, avemag_bin[i], m1_bin[i], m2_bin[i], sigma_e1_bin[i], sigma_e2_bin[i], ngal_bin[i]

