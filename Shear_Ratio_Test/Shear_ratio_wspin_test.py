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
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

DFLAG ='' # !!! TEMPORARY FLAG TO CALC DATA GT FROM MICE OCTANT, NOT FIDUCIAL MICE.

Cov_Method = "Spin" #"Spin"  # The method for calculating the gamma_t realisations for use in covariance estimation
                     # "Spin" means do many spin realisations of the source ellipticities (ie - shape noise only)
                     # "Patch" means read in the other MICE realisations (1/8 of the sky)
                     # divide them into patches and calcute the gamma_t from each patch.
if Cov_Method != "Spin" and Cov_Method != "Patch" and Cov_Method != "None":
    print("Cov_Method must be set to Spin, Patch or None. Currently it's set to %s. EXITING." %Cov_Method)
    sys.exit()
    
nPatch = 16           # If Cov_Method is Patch, the MICE octant is split into nPatch RA
                     # and nPatch Dec slices (nPatch^2 patches in total). gamma_t is
                     # calculated from each patch.

nlens = 5            # !!! IMPORTANT PARAMETER.
                     # The number of bins the lenses are divided into.
                     
                     
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
    paramfile = sys.argv[8]

# The following lines are used to separately save a finely binned gamma_t
# for use in producing theoretical gamma_t's at the correct (mean-point) theta
# values using kcap.
if nbins > 10:
    thsavetag = '_FineBins%s' %nbins
else:
    thsavetag =	''
    

# ----- Load Input ----- #
from Get_Input import Get_Input              
GI = Get_Input(paramfile)
SOURCE_TYPE = GI.Source_Type()
LENS_TYPE = GI.Lens_Type()     
RANDOM_TYPE = GI.Random_Type() 

if SOURCE_TYPE == "K1000":    
    Blind = GI.Blind()
    CATDIR='/disk09/KIDS/K1000_TWO_PT_STATS/TOMOCATS/'

    # These are the old catalogues binned using DIR-redshift estimates
    #fitscat=CATDIR+'/K1000_N_BLIND_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_goldclasses_5Z_%s.fits' %(Blind,JZBIN)

    # These are the new catalogues binned using the SOM-redshift estimates
    SOMFLAGNAME = GI.SOMFLAGNAME()
    fitscat=CATDIR+'K1000_N_BLIND_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_%s_5Z_%s.fits' %(Blind,SOMFLAGNAME,JZBIN)
    
    # Fits extension and keywords depend on data source:
    iext=1     # File extension
    ra_keyword='ALPHA_J2000'
    dec_keyword='DELTA_J2000'
    w_keyword='weight'
    e1_keyword='e1'
    e2_keyword='e2'
    flip_g1 = False
    flip_g2 = False

    # Extra info needed to identify lenscat
    lens_tag = ''
    # Location of the output files 
    OUTDIR='Output/SOURCE-%s_LENS-%s_Blind%s_SOM%s' %(SOURCE_TYPE,LENS_TYPE,Blind,SOMFLAGNAME)
    
elif SOURCE_TYPE == "MICE2_KV450":
    Mag_OnOff = GI.Mag_OnOff()
    Pz = GI.Pz_TrueEstimated()
    SN = GI.SN()
    MICE_DIR = '/home/bengib/MICE2_Mocks/MICE2_KV450/shear_CATS_Pz%s_SN%s%s' %(Pz,SN,DFLAG)
    fitscat='%s/MICE2_KV450_magnification_%s_small_6Z_%s.fits' %(MICE_DIR,Mag_OnOff,JZBIN)

    iext=1
    ra_keyword='ra_gal'
    dec_keyword='dec_gal'
    w_keyword='recal_weight'
    e1_keyword='gamma1'
    e2_keyword='gamma2'
    if Mag_OnOff == 'on':
        ra_keyword += '_mag'
        dec_keyword += '_mag'
    flip_g1 = True   # Quite bizarrely one must flip g1 with Jan Luca's MICE mocks
    flip_g2 = False  # (But not with Chris Blake's MICE catalogues. Simples.)

    # Extra info needed to identify lenscat
    lens_tag = '_mag%s'%Mag_OnOff
    # Location of the output files
    OUTDIR='Output/SOURCE-%s_LENS-%s_Pz%s_SN%s_mag%s' %(SOURCE_TYPE,LENS_TYPE,Pz,SN,Mag_OnOff)
    
else:
    print("This code only accepts Source_Type set to 'K1000' or 'MICE2_KV450'... ")
    print("Here SOURCE_TYPE is set to ",SOURCE_TYPE," which is not recognised so I'm exiting.")
    sys.exit(1)
    

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
    os.makedirs(OUTDIR+'/SPIN')

# Because we're rotating the observed ellipticities of the sources
# for the spin test we won't use the built
# in treecorr read from the fits catalogue, but read it first and then modify

#open the fits catalogue
fitsfile = fits.open(fitscat)

# read in position and shear
ra = (fitsfile[iext].data[ra_keyword])
dec = (fitsfile[iext].data[dec_keyword])
eobs1 = (fitsfile[iext].data[e1_keyword])
eobs2 = (fitsfile[iext].data[e2_keyword])
weight = (fitsfile[iext].data[w_keyword])

ngals=len(ra)
print("ngals", ngals)

# the unspun source catalogue
sourcecat = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1,g2=eobs2,
                             ra_units='deg', dec_units='deg',w=weight,
                             flip_g1=flip_g1, flip_g2=flip_g2)
sourcecat_wsq = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1,g2=eobs2,
                             ra_units='deg', dec_units='deg',w=weight**2,
                             flip_g1=flip_g1, flip_g2=flip_g2)

# If we're using the MICE octant for covariance, read in the source info here
if Cov_Method == "Patch":
    MICE_DIR_OCTANT = '/home/bengib/MICE2_Mocks/MICE2_KV450/shear_CATS_Pz%s_SN%s_Octant' %(Pz,SN)
    fitsfile_o=fits.open('%s/MICE2_KV450_magnification_%s_small_6Z_%s.fits' %(MICE_DIR_OCTANT,Mag_OnOff,JZBIN))
    # read in position and shear for octant sources
    ra_s = (fitsfile_o[iext].data[ra_keyword])
    dec_s = (fitsfile_o[iext].data[dec_keyword])
    eobs1_s = (fitsfile_o[iext].data[e1_keyword])
    eobs2_s = (fitsfile_o[iext].data[e2_keyword])
    weight_s = (fitsfile_o[iext].data[w_keyword])

# Useful functions if Cov_Method is Patch
# Function returns elements INSIDE (ra,dec) patch
def Select_Patch(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the elements in Q corresponding to INSIDE the (ra,dec) range 
    idx_ra = np.where(np.logical_and(ra<rhi, ra>rlo))[0]
    idx_dec = np.where(np.logical_and(dec<dhi, dec>dlo))[0]
    idx = np.intersect1d(idx_ra, idx_dec) 
    return Q[idx]
# Function returns elements OUTSIDE (ra,dec) patch
def Split_Fields(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the elements in Q corresponding to OUTSIDE the (ra,dec) range (Jackknife)
    idx_ra = np.append( np.where(ra>rhi)[0], np.where(ra<rlo)[0] )
    idx_dec = np.append( np.where(dec>dhi)[0], np.where(dec<dlo)[0] )
    idx = np.unique( np.append(idx_ra, idx_dec) ) # delete the things that appear more than once
    return Q[idx]




# To minimise I/O source catalogue read across the cluster, we'll loop over all lens bins
# If the WL-72 nodes workers are not available this will be sub-optimal and it would be better
# to queue each source-lens bin pair individually on the cluster

for IZBIN in range (izin,izin+1):   #1,2,3,4,5

    lenscatname='LENSCATS/%s%s%s/lens_cat_%sZ_%s.fits' %(LENS_TYPE,lens_tag,DFLAG, nlens,IZBIN)
    rancatname='LENSCATS/%s%s%s/lens_cat_%sZ_%s.fits' %(RANDOM_TYPE,lens_tag,DFLAG, nlens,IZBIN)
    outfile_main='%s/GT_6Z_source_%s_%sZ_lens_%s%s%s.asc' %(OUTDIR,JZBIN, nlens,IZBIN, DFLAG,thsavetag)

    # the lens catalogue we will not modify so we can use the built in treecorr option to read 
    # in directly from the catalogue
    lenscat = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg',
                               w_col='WEICOMP')
    rancat = treecorr.Catalog(rancatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg',
                              w_col='WEICOMP')


    # Set-up the different correlations that we want to measure
    bin_slop=0.08 # optimised in Flinc sims
    #bin_slop=0.12 # faster option

    # Number of source lens pairs
    nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',
                                  bin_slop=bin_slop)
    # Number of source random pairs     
    nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',
                                  bin_slop=bin_slop)
    # Average shear around lenses     
    ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',
                                bin_slop=bin_slop)
    # Average shear around randoms     
    rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',
                                bin_slop=bin_slop)

    # ----- Average around lenses with squared weights -----
    f = fits.open(lenscatname)
    ra_l_wsq = f[1].data['ALPHA_J2000']
    dec_l_wsq = f[1].data['DELTA_J2000']
    w_l_wsq = f[1].data['WEICOMP']
    f.close()
    lenscat_wsq = treecorr.Catalog(ra=ra_l_wsq,dec=dec_l_wsq,
                                   ra_units='deg', dec_units='deg',
                                   w=w_l_wsq**2 )
    ls_wsq = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin',
                                    bin_slop=bin_slop)
    # -------------------------------------------------------

    # Now calculate the different 2pt correlation functions
    print("Calculating the number of source-lens pairs")
    nlns.process(lenscat,sourcecat)
    print("Calculating the number of source-random pairs")
    nrns.process(rancat,sourcecat)
    print("Calculating the average shear around the lenses")
    ls.process(lenscat,sourcecat)    # only this one needs to be recalculated for each spin test
    print("Calculating the average shear around the randoms")
    rs.process(rancat,sourcecat)
    print("Calculating the average shear around the lenses WITH SQUARED WEIGHTS")
    ls_wsq.process(lenscat_wsq,sourcecat_wsq)


    

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
    # This next line produces data used in the computation of the analytical covariance matrix
    #npairs_weighted_simple = (ls.weight)*(ls.weight)/(ls_wsq.weight)
    print("Writing out", outfile_main)
    
    #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
    treecorr.util.gen_write(outfile_main,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
             'rangamT','rangamX','ransigma', 'ranweight','weight_sqrd' ],
            [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), ls.weight, ls.npairs,
              ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi), rs.weight,ls_wsq.weight ])

    if Cov_Method == "Spin":
        # WE ARE CALCULATING THE COVARIANCE BY
        # spinning the source ellipticities many times (shape noise only)
    
        # Now carry out the spin test with default binslop
        # Average shear around lenses with default fast binslop    
        lssp = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin')
        # ----- 19/12/19: GIBLIN's EDIT - CALC RANDOM SIGNAL AS WELL FOR SOURCE-SPINS ----- 
        rssp = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin')
        # ----- #
    
        # As we loop through lens bins, we want to have the same random spins of the sources
        # for each lens bin. Therefore we have to fix the seed.
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

            sourcespin = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1_spin,g2=eobs2_spin,
                                          ra_units='deg', dec_units='deg',w=weight,
                                          flip_g1=flip_g1, flip_g2=flip_g2)
            lssp.process(lenscat,sourcespin)    # only this needs to be recalculated for each spin test

            # ----- 19/12/19: GIBLIN's EDIT - CALC RANDOM SIGNAL AS WELL FOR SOURCE-SPINS -----
            # ----- (MIGHT BE NECESSARY FOR MICE MOCKS) -----
            rssp.process(rancat,sourcespin)
            # ----- #

            # lets write these all out and post-process because we don't know how many we will need in the end
            # note it another set of trials is run - change the fixed seed number

            #in this case however we do not subtract the random signal as the shears are randomised
            #and therefore the random signal should be zero within the noise
            #in principle rs should also be calculated but we do not include this extra noise term
            #to speed up the spin test calculation
            gamma_t = lssp.xi*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rssp.xi         # GIBLIN! NOTE YOU CHANGED THIS
            gamma_x = lssp.xi_im*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rssp.xi_im   # rs --> rssp!!!

            outfile_spin='%s/SPIN/GT_SPIN_%s_6Z_source_%s_%sZ_lens_%s%s.asc' %(OUTDIR,i,JZBIN, nlens,IZBIN,thsavetag)

            #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
            treecorr.util.gen_write(outfile_spin,
                            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
                             'rangamT','rangamX','ransigma' ],
                            [ lssp.rnom,lssp.meanr, lssp.meanlogr,gamma_t, gamma_x,
                              np.sqrt(lssp.varxi), lssp.weight, lssp.npairs,
                              lssp.xi, lssp.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])

        
            # We have a memory leak when doing these spin realisations (memory build up until job dies).
            # Hopefully this line will stop this from happening:
            sourcespin.gfields.clear()
        
            # Let's attempt to stop this by deleting all of these objects at the end of each spin:
            #del theta, ct, st, eobs1_spin, eobs2_spin
            #del sourcespin, gamma_t, gamma_x, outfile_spin


    elif Cov_Method == "Patch":
        # WE ARE CALCULATING THE COVARIANCE BY
        # Reading in the MICE octant, dividing it up and calculating gamma_t
        # from the segments.

        # Octant sources have already been read in. Just need to grab the lens data.
        # This needs to be read in so we can divide the octant up.
        lensfile_o = fits.open( 'LENSCATS/%s%s_Octant/lens_cat_%sZ_%s.fits' %(LENS_TYPE,lens_tag, nlens,IZBIN) )
        ranfile_o = fits.open( 'LENSCATS/%s%s_Octant/lens_cat_%sZ_%s.fits' %(RANDOM_TYPE,lens_tag, nlens,IZBIN) )
                
        # Lenses
        ra_l =  (lensfile_o[iext].data['ALPHA_J2000'])
        dec_l = (lensfile_o[iext].data['DELTA_J2000'])
        w_l =   (lensfile_o[iext].data['WEICOMP'])
        # Randoms
        ra_r =  (ranfile_o[iext].data['ALPHA_J2000'])
        dec_r = (ranfile_o[iext].data['DELTA_J2000'])
        w_r =   (ranfile_o[iext].data['WEICOMP'])
        # (Sources read in prior to loop)

        # Boundaries of the RA,Dec patches
        ra_coarse = np.linspace(0., 90., nPatch+1)
        dec_coarse = np.linspace(0., 90., nPatch+1)

        # The functions for recalculating the lens/random gamma_t's in each patch
        bin_slop_fast = 0.12
        lssp = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=bin_slop_fast)
        rssp = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=bin_slop_fast)

        # Cycle through the Patches
        for i in range(int(nPatch*nPatch)):
            
            print("On lens %s, patch %s of %s" %(IZBIN,i,nPatch*nPatch-1))
            ridx = i % nPatch
            didx = int( i / nPatch )
            rlo = ra_coarse[ridx]
            rhi = ra_coarse[ridx+1]
            dlo = dec_coarse[didx]
            dhi = dec_coarse[didx+1]
            # Find sources in the patch
            tmp_ra_s = Split_Fields(ra_s,    ra_s, dec_s, rlo, rhi, dlo, dhi)
            tmp_dec_s = Split_Fields(dec_s,  ra_s, dec_s, rlo, rhi, dlo, dhi)
            tmp_e1_s = Split_Fields(eobs1_s, ra_s, dec_s, rlo, rhi, dlo, dhi)
            tmp_e2_s = Split_Fields(eobs2_s, ra_s, dec_s, rlo, rhi, dlo, dhi)
            tmp_w_s = Split_Fields(weight_s, ra_s, dec_s, rlo, rhi, dlo, dhi)
            # Find lenses in the patch
            tmp_ra_l = Split_Fields(ra_l,    ra_l, dec_l, rlo, rhi, dlo, dhi)
            tmp_dec_l = Split_Fields(dec_l,  ra_l, dec_l, rlo, rhi, dlo, dhi)
            tmp_w_l = Split_Fields(w_l,      ra_l, dec_l, rlo, rhi, dlo, dhi)
            # Find randoms in the patch
            tmp_ra_r = Split_Fields(ra_r,    ra_r, dec_r, rlo, rhi, dlo, dhi)
            tmp_dec_r = Split_Fields(dec_r,  ra_r, dec_r, rlo, rhi, dlo, dhi)
            tmp_w_r = Split_Fields(w_r,      ra_r, dec_r, rlo, rhi, dlo, dhi)

            # Set up the catalogues:
            tmp_sourcecat = treecorr.Catalog(ra=tmp_ra_s, dec=tmp_dec_s,
                                             g1=tmp_e1_s, g2=tmp_e2_s,
                                             ra_units='deg', dec_units='deg',w=tmp_w_s,
                                             flip_g1=True, flip_g2=False)
            tmp_lenscat = treecorr.Catalog(ra=tmp_ra_l, dec=tmp_dec_l,
                                             ra_units='deg', dec_units='deg',w=tmp_w_l)
            tmp_rancat = treecorr.Catalog(ra=tmp_ra_r, dec=tmp_dec_r,
                                             ra_units='deg', dec_units='deg',w=tmp_w_r)

            lssp.process(tmp_lenscat,tmp_sourcecat)    
            rssp.process(tmp_rancat,tmp_sourcecat)

            # lets write these all out and post-process because we don't know how many we will need in the end
            # Here we do subtract the random signal 
            gamma_t = lssp.xi*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rssp.xi
            gamma_x = lssp.xi_im*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rssp.xi_im

            if not os.path.exists(OUTDIR+'/PATCH/'):
                os.makedirs(OUTDIR+'/PATCH/')
            tmp_outfile='%s/PATCH/GT_PATCH_%sof%s_6Z_source_%s_%sZ_lens_%s%s.asc' %(OUTDIR,i,int(nPatch*nPatch),JZBIN,nlens,IZBIN, thsavetag)

            #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
            treecorr.util.gen_write(tmp_outfile,
                            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX',
                             'rangamT','rangamX','ransigma' ],
                            [ lssp.rnom,lssp.meanr, lssp.meanlogr,gamma_t, gamma_x,
                              np.sqrt(lssp.varxi), lssp.weight, lssp.npairs,
                              lssp.xi, lssp.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ])


            # Lines to avoid memory building up and crashing:
            tmp_sourcecat.gfields.clear()
            tmp_lenscat.nfields.clear()
            tmp_rancat.nfields.clear()
