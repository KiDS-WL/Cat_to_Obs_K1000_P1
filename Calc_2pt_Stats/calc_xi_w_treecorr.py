# ----------------------------------------------------------------
# File Name:           calc_xi_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate xi_+/- 
#                      given a KiDS fits catalogue
#                      script will need to change if keywords in KIDS cats are updated
# ----------------------------------------------------------------

import treecorr
import sys


# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, fitscat1, fitscat2, outfilename
if len(sys.argv) <7: 
    #print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log? catalogue1.fits \
    #        catalogue2.fits outfilename" % sys.argv[0]) 
    #print "Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? \
    #       catalogue1.fits catalogue2.fits outfilename" % sys.argv[0] 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 
    lin_not_log = sys.argv[4] 
    fitscat1 = sys.argv[5]
    fitscat2 = sys.argv[6]
    outfile = sys.argv[7]

# prepare the catalogues

cat1 = treecorr.Catalog(fitscat1, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')
cat2 = treecorr.Catalog(fitscat2, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')

# Define the binning based on command line input
# Using bin_slop of 0.08 which Linc found to be optimal for xi_+/- on the mocks
if(lin_not_log=='true'): 
    gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
         bin_type='Linear', bin_slop=0.08)
else:
    #gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
    #bin_type='Log', bin_slop=0.08)
    gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=0.08)

# Calculate the 2pt correlation function
gg.process(cat1,cat2)

# Write it out to a file and praise-be for Jarvis and his well documented code
gg.write(outfile)