# ----------------------------------------------------------------
# File Name:           calc_gt_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate gamma_t
#                      we're using the Mandelbaum estimator where the randoms are subtracted
#                      so we need a lens, source and random catalogue 
#                      script will need to change if keywords in KIDS cats are updated
#  Treecorr doc: https://rmjarvis.github.io/TreeCorr/_build/html/correlation2.html
# ----------------------------------------------------------------

import treecorr
import sys
import numpy as np


# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, lenscat, rancat, sourcecat, outfilename
if len(sys.argv) <8: 
    #print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? lenscat.fits \
    #        randomcat.fits sourcecat.fits outfilename" % sys.argv[0]) 
    print "Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? \
           lenscat.fits randomcat.fits sourcecat.fits catalogue1.fits catalogue2.fits outfilename" % sys.argv[0] 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 
    lin_not_log = sys.argv[4] 
    lenscatname = sys.argv[5]
    rancatname = sys.argv[6]
    sourcecatname = sys.argv[7]
    outfile = sys.argv[8]

# prepare the catalogues

lenscat = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')
rancat = treecorr.Catalog(rancatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')
sourcecat = treecorr.Catalog(sourcecatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')

# Define the binning based on command line input
# Using bin_slop of 0.08 which Linc found to be optimal for xi_+/- on the mocks
if(lin_not_log=='true'): 
    gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
         bin_type='Linear', bin_slop=0.08)
else:
    #gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
    #bin_type='Log', bin_slop=0.08)

# Set-up the different correlations that we want to measure

    # Number of source lens pairs
    nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=0.08)
    # Number of source random pairs     
    nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=0.08)
    # Average shear around lenses     
    ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=0.08)
    # Average shear around randoms     
    rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=0.08)

# Now calculate the different 2pt correlation functions
nlns.process(lenscat,sourcecat)
nrns.process(rancat,sourcecat)
ls.process(lenscat,sourcecat)
rs.process(rancat,sourcecat)

# We will use the Mandelbaum 2006 estimator which includes both the random and boost correction.
# It is given by
# gt = (SD-SR)/NsNr
# SD = sum of shear around source-lens pairs
# SR = sum of shear around source-random pairs
# NrNs = number of source random pairs
# Note that we have a different number of randoms from lenses so we also need to make
# sure that we rescale NrNs accordingly

# The easiest way to do this in Treecorr is
# gt = (SD/NlNs)*NlNs/NrNs - SR/NrNs
# where NsNl = number of source lens pairs
# 
# SD/NsNl = ls.xi
# NlNs = nlns.weight/nlns.tot
# NrNs = nrns.weight/nrns.tot
# SR/NrNs = rs.xi

gamma_t = ls.xi*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi
gamma_x = ls.xi_im*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi_im

#Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
treecorr.util.gen_write(outfile,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
            'rangamT','rangamX','ransigma' ],
            [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), ls.weight, ls.npairs,
            ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])