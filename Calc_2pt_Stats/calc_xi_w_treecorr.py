# ----------------------------------------------------------------
# File Name:           calc_xi_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate xi_+/- 
#                      given a KiDS fits catalogue
#                      script will need to change if keywords in KIDS cats are updated
# ----------------------------------------------------------------

import treecorr
import sys
import numpy as np

# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, fitscat1, fitscat2, outfilename
if len(sys.argv) <9: 
    print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log? catalogue1.fits \
            catalogue2.fits outfilename weighted_analysis?" % sys.argv[0]) 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 
    lin_not_log = sys.argv[4] 
    fitscat1 = sys.argv[5]
    fitscat2 = sys.argv[6]
    outfile = sys.argv[7]
    weighted = sys.argv[8]
    
if len(sys.argv) >9:  # optional column names supplied
    cat1e1name = sys.argv[9]
    cat1e2name = sys.argv[10]
    cat2e1name = sys.argv[11]
    cat2e2name = sys.argv[12]
    fine_binning = sys.argv[13]
else: #use defaults
    cat1e1name = 'e1'
    cat1e2name = 'e2'
    cat2e1name = 'e1'
    cat2e2name = 'e2'
    fine_binning = 'true'

if fine_binning:
    # when using fine bins I find this is suitable
    inbinslop = 1.5
else:
    # when using broad bins it needs to be much finer
    inbinslop = 0.08

# prepare the catalogues

cat1 = treecorr.Catalog(fitscat1, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col=cat1e1name, g2_col=cat1e2name, w_col='weight')
cat2 = treecorr.Catalog(fitscat2, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col=cat2e1name, g2_col=cat2e2name, w_col='weight')

# Define the binning based on command line input
if(lin_not_log=='true'): 
    gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
         bin_type='Linear', bin_slop=inbinslop)
else: # Log is the default bin_type for Treecorr
    gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop)

# Calculate the 2pt correlation function
gg.process(cat1,cat2)

if (weighted=='true'):    
# prepare the weighted_square catalogues - hack so that Treecorr returns the correct Npairs for a weighted sample

    cat1_wsq = treecorr.Catalog(fitscat1, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col=cat1e1name, g2_col=cat1e2name, w_col='weightsq')
    cat2_wsq = treecorr.Catalog(fitscat2, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col=cat2e1name, g2_col=cat2e2name, w_col='weightsq')

    # Define the binning based on command line input
    if(lin_not_log=='true'): 
        gg_wsq = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',\
                                    bin_type='Linear', bin_slop=inbinslop)
    else: # Log is the default bin_type for Treecorr
        gg_wsq = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                    bin_slop=inbinslop)    

    # Calculate the weighted square 2pt correlation function
    gg_wsq.process(cat1_wsq,cat2_wsq)

    # Calculate the weighted Npairs = sum(weight_a*weight_b)^2 / sum(weight_a^2*weight_b^2)
    
    npairs_weighted = (gg.weight)*(gg.weight)/gg_wsq.weight

    #Use treecorr to write out the output file updating the npairs column and praise-be for Jarvis and his well documented code                                                   
    treecorr.util.gen_write(outfile,
            ['r_nom','meanr','meanlogr','xip','xim','xip_pm','xim_im','sigma_xip', 'npairs', 'weight','npairs_weighted' ],
            [ gg.rnom,gg.meanr, gg.meanlogr,gg.xip, gg.xim, gg.xip_im, gg.xim_im, np.sqrt(gg.varxip), gg.npairs, 
            gg.weight, npairs_weighted])                                                  
else:    
    # Write it out unweighted npairs and praise-be again for Jarvis and his well documented code
    gg.write(outfile)
