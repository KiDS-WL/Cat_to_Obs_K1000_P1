# ----------------------------------------------------------------
# File Name:           Shear_ratio_wspin_test.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate GGL 
#                      for the shear ratio test
#                      for the covariance we use the spin test
#                      where the observed ellipticity is randomised
# ----------------------------------------------------------------

import treecorr
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, fitscat1, fitscat2, outfilename
if len(sys.argv) <6: 
    #print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) source_tomobin number_of_spins" % sys.argv[0]) 
    print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) source_tomobin" % sys.argv[0]) 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 
    JZBIN = int(sys.argv[4]) 
    ntrials = int(sys.argv[5]) 
    izin = int(sys.argv[6]) 
    ntprevrun = int(sys.argv[7]) 

Blind='A'
# Source Catalogues
CATDIR='/disk09/KIDS/K1000_TWO_PT_STATS/'
fitscat=CATDIR+'/TOMOCATS/K1000_N_BLIND_'+Blind+'_v3_6Z_'+str(JZBIN)+'.fits'

# Location of the output files
OUTDIR='/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/Shear_Ratio_Test/Output/'
#CATDIR+'/OUTSTATS/SHEAR_RATIO/GT/'

# Because we're rotating the observed ellipticities of the sources
# for the spin test we won't use the built
# in treecorr read from the fits catalogue, but read it first and then modify

#open the fits catalogue
fitsfile = fits.open(fitscat)

# read in position and shear
ra = (fitsfile[1].data['ALPHA_J2000'])
dec = (fitsfile[1].data['DELTA_J2000'])
eobs1 = (fitsfile[1].data['e1'])
eobs2 = (fitsfile[1].data['e2'])
weight = (fitsfile[1].data['weight'])

ngals=len(ra)
print("ngals", ngals)

# the unspun source catalogue
sourcecat = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1,g2=eobs2,ra_units='deg', dec_units='deg',w=weight)

# To minimise I/O source catalogue read across the cluster, we'll loop over all lens bins
# If the WL-72 nodes workers are not available this will be sub-optimal and it would be better
# to queue each source-lens bin pair individually on the cluster

for IZBIN in range (izin,izin+1):   #1,2,3,4,5

    lenscatname=CATDIR+'/GGLCATS/BOSS_data_5Z_'+str(IZBIN)+'.fits'
    rancatname=CATDIR+'/GGLCATS/BOSS_random_5Z_'+str(IZBIN)+'.fits'
    outfile_main=OUTDIR+'/K1000_GT_6Z_source_'+str(JZBIN)+'_5Z_lens_'+str(IZBIN)+'.asc'

    # the lens catalogue we will not modify so we can use the built in treecorr option to read 
    # in directly from the catalogue
    lenscat = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')
    rancat = treecorr.Catalog(rancatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  w_col='WEICOMP')


    # Set-up the different correlations that we want to measure
    bin_slop=0.08 # optimised in Flinc sims
    #bin_slop=0.12 # faster option

    # Number of source lens pairs
    nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
        bin_slop=bin_slop)
    # Number of source random pairs     
    nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
        bin_slop=bin_slop)
    # Average shear around lenses     
    ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
        bin_slop=bin_slop)
    # Average shear around randoms     
    rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
        bin_slop=bin_slop)

    # Now calculate the different 2pt correlation functions
    print("Calculating the number of source-lens pairs")
    nlns.process(lenscat,sourcecat)
    print("Calculating the number of source-random pairs")
    nrns.process(rancat,sourcecat)
    print("Calculating the average shear around the lenses")
    ls.process(lenscat,sourcecat)    # only this one needs to be recalculated for each spin test
    print("Calculating the average shear around the randoms")
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

    print("Writing out", outfile_main)
    
    #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
    treecorr.util.gen_write(outfile_main,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
            'rangamT','rangamX','ransigma' ],
            [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), ls.weight, ls.npairs,
            ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])


    # Now carry out the spin test with default binslop
    # Average shear around lenses with default fast binslop    
    lssp = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin')

    # As we loop over lens bins and we want to make sure that we have the same spin
    # for each of the source trials analysis of each lens bin we have to fix the seed
    np.random.seed(42)

    print("Running spin test using ntrials =", ntrials)

    if ntprevrun > 0 :
        for i in range (ntprevrun):
            # spin the last galaxy sample
            theta = np.random.uniform(0,np.pi,ngals)

    for i in range (ntprevrun,ntrials):
        # spin the last galaxy sample
        theta = np.random.uniform(0,np.pi,ngals)
        ct= np.cos(2.0*theta)
        st= np.sin(2.0*theta)
        eobs1_spin=eobs1*ct + eobs2*st
        eobs2_spin=-1.0*eobs1*st + eobs2*ct

        sourcespin = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1_spin,g2=eobs2_spin,ra_units='deg', dec_units='deg',w=weight)
        lssp.process(lenscat,sourcespin)    # only this needs to be recalculated for each spin test

        # lets write these all out and post-process because we don't know how many we will need in the end
        # note it another set of trials is run - change the fixed seed number

        #in this case however we do not subtract the random signal as the shears are randomised
        #and therefore the random signal should be zero within the noise
        #in principle rs should also be calculated but we do not include this extra noise term
        #to speed up the spin test calculation
        gamma_t = lssp.xi*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi
        gamma_x = lssp.xi_im*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi_im


        outfile_spin=OUTDIR+'/SPIN/K1000_GT_SPIN_'+str(i)+'_6Z_source_'+str(JZBIN)+'_5Z_lens_'+str(IZBIN)+'.asc'

        #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
        treecorr.util.gen_write(outfile_spin,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
            'rangamT','rangamX','ransigma' ],
            [ lssp.rnom,lssp.meanr, lssp.meanlogr,gamma_t, gamma_x, np.sqrt(lssp.varxi), lssp.weight, lssp.npairs,
            lssp.xi, lssp.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])

