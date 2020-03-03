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
if len(sys.argv) < 10: 
    print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? \
          lenscat.fits randomcat.fits sourcecat.fits catalogue1.fits catalogue2.fits outfilename" % sys.argv[0])
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
    weighted = sys.argv[9]

# prepare the catalogues

lenscat = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')
rancat = treecorr.Catalog(rancatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')
sourcecat = treecorr.Catalog(sourcecatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')


#to do include wcs cut
#,flag_col='KIDSMASK',ignore_flag=16384

# the value of bin_slop
fine_binning = True

if fine_binning:
    # when using fine bins I find this is suitable
    inbinslop_NN = 1.0
    inbinslop_NG = 1.2
else:
    # when using broad bins it needs to be much finer
    inbinslop_NN = 0.03
    inbinslop_NG = 0.05


# Define the binning based on command line input
if(lin_not_log=='true'): 
    # Number of source lens pairs
    nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_type='Linear', bin_slop=inbinslop_NN)
    # Number of source random pairs     
    nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_type='Linear', bin_slop=inbinslop_NN)
    # Average shear around lenses     
    ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_type='Linear', bin_slop=inbinslop_NG)
    # Average shear around randoms     
    rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_type='Linear', bin_slop=inbinslop_NG)

else:
    # Set-up the different correlations that we want to measure

    # Number of source lens pairs
    nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop_NN)
    # Number of source random pairs     
    nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop_NN)
    # Average shear around lenses     
    ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop_NG)
    # Average shear around randoms     
    rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop_NG)

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

if (weighted=='true'):   
    # prepare the weighted_square catalogues - hack so that Treecorr returns the correct Npairs for a weighted sample
    
    lenscat_wsq = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                   w_col='WEICOMPsq')
    sourcecat_wsq = treecorr.Catalog(sourcecatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                   g1_col='e1', g2_col='e2', w_col='weightsq')

    # Define the binning based on command line input
    if(lin_not_log=='true'): 
        # Average shear around lenses     
        ls_wsq = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                        bin_type='Linear', bin_slop=inbinslop_NG)
    else: # Log is the default bin_type for Treecorr   
        ls_wsq = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                        bin_slop=inbinslop_NG)

    # Calculate the weighted square 2pt correlation function
    ls_wsq.process(lenscat_wsq,sourcecat_wsq)
    # Calculate the weighted Npairs = sum(weight_a*weight_b)^2 / sum(weight_a^2*weight_b^2)
    
    npairs_weighted = (ls.weight)*(ls.weight)/ls_wsq.weight

    #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
    treecorr.util.gen_write(outfile,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs_weighted', 'nocor_gamT', 'nocor_gamX', 
            'rangamT','rangamX','ransigma' ],
            [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), npairs_weighted, ls.npairs,
            ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])

else:     
     #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
     treecorr.util.gen_write(outfile,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
            'rangamT','rangamX','ransigma' ],
            [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), ls.weight, ls.npairs,
            ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])
