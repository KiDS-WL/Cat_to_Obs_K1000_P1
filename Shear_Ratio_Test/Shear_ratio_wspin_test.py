# ----------------------------------------------------------------
# File Name:           metacal_mock_calc_xi_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate xi_+/- 
#                      given a KiDS mock fits catalogue
#                      values of m are added and the intrinsic ellipticity is randomisd
#                      so we can test neff with metacal weights
# ----------------------------------------------------------------

import treecorr
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, fitscat1, fitscat2, outfilename
if len(sys.argv) <3: 
    #print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log? catalogue1.fits \
    #        catalogue2.fits outfilename" % sys.argv[0]) 
    #print "Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? \
    #       catalogue1.fits catalogue2.fits outfilename" % sys.argv[0] 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 

fitscat="/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/G9_N_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321.cat"

# Because we're adding m's and rotating the intrinsic ellipticities we won't use the built
# in treecorr read from the fits catalogue, but read it first and then modify
#open the fits catalogue
fitsfile = fits.open(fitscat)
area=54*60*60.0  # rough G9 field size

# read in position and shear
ra = (fitsfile[2].data['ALPHA_J2000'])
dec = (fitsfile[2].data['DELTA_J2000'])
eobs1 = (fitsfile[2].data['raw_e1'])
eobs2 = (fitsfile[2].data['raw_e2'])
Rm = 1+(fitsfile[2].data['metacal_m1'])

#fix the metacal m and weight - currently random between 0 and 1/10
ngals=len(ra)
print "ngals", ngals
w = np.ones(ngals)
weight_w_m = w*Rm

# m-calibration corrections for two alternative weights
Kcor=(np.sum(w*Rm)/np.sum(w))**2
Kcor_wm=(np.sum(weight_w_m*Rm)/np.sum(weight_w_m))**2

print "Kcor w=1: %4.2f Kcor w=(1+m): %4.2f" %(Kcor, Kcor_wm)

# calculate neffective
neff_orig =((np.sum(w))**2/np.sum(w*w))/area
neff_new = ((np.sum(w*Rm))**2/np.sum(w*Rm*w*Rm))/area
neff_new_wm = ((np.sum(weight_w_m*Rm))**2/np.sum(weight_w_m*Rm*weight_w_m*Rm))/area

print "neff orig: %4.2f neff new: %4.2f neff new_wm: %4.2f" %(neff_orig, neff_new, neff_new_wm)

# and sigma_e
sige_orig = np.sum(w*w*(eobs1*eobs1 + eobs2*eobs2))/np.sum(w*w)
sige_new = np.sum(w*w*(eobs1*eobs1 + eobs2*eobs2))/np.sum(w*Rm*w*Rm)
sige_new_wm = np.sum(weight_w_m*weight_w_m*(eobs1*eobs1 + eobs2*eobs2))/np.sum(weight_w_m*Rm*weight_w_m*Rm)

print "sigma_e: orig: %4.2f new: %4.2f new_wm: %4.2f" % (np.sqrt(sige_orig), np.sqrt(sige_new),np.sqrt(sige_new_wm))

# Now lets do the spin test 
ntrials = 30
xipspin=np.zeros(ntrials)
xipspin_wm=np.zeros(ntrials)
inbinslop=0.1

for i in range (ntrials):
    # spin the last galaxy sample
    theta = np.random.uniform(0,np.pi,ngals)
    ct= np.cos(2.0*theta)
    st= np.sin(2.0*theta)
    eobs1_spin=eobs1*ct + eobs2*st
    eobs2_spin=-1.0*eobs1*st + eobs2*ct

    catspin = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1_spin,g2=eobs2_spin,ra_units='deg', dec_units='deg',w=w)
    catspin_wm = treecorr.Catalog(ra=ra,dec=dec,g1=eobs1_spin,g2=eobs2_spin,ra_units='deg', dec_units='deg',w=weight_w_m)

    ggspin = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop)
    ggspin_wm = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
         bin_slop=inbinslop)

    ggspin.process(catspin)
    ggspin_wm.process(catspin_wm)
    xipspin[i]=ggspin.xip[0]/Kcor
    xipspin_wm[i]=ggspin_wm.xip[0]/Kcor_wm

print "mean spin:%.3e spin_wm:%.3e" % (np.average(xipspin),np.average(xipspin_wm))
print "std  spin: %.3e spin_wm: %.3e" % (np.std(xipspin),np.std(xipspin_wm))

#Npairs=ggmw.npairs  # does not include weights - this is a factor of 2 smaller than the Npair_theory definition
theta=ggspin.meanr
#Area = ngals*ngals*2*np.pi*theta*(theta_max-theta_min)/(Npairs*2)
#print Area/(60*60)  # yup recovers expected area

# predicted noise covariance
Cov_orig = sige_orig**2 / (area*neff_orig**2*2*np.pi*theta[0]*(theta_max-theta_min))
Cov_new = sige_new**2 / (area*neff_new**2*2*np.pi*theta[0]*(theta_max-theta_min))
Cov_new_wm = sige_new_wm**2 / (area*neff_new_wm**2*2*np.pi*theta[0]*(theta_max-theta_min))
Cov_mix = sige_new**2 / (area*neff_orig**2*2*np.pi*theta[0]*(theta_max-theta_min))

print "Cov: new: %.3e new_wm: %.3e orig: %.3e mix: %.3e " % (np.sqrt(Cov_new),np.sqrt(Cov_new_wm),np.sqrt(Cov_orig),np.sqrt(Cov_mix))

#    print ggex.xip, ggobs.xip/Kcor,ggmw.xip/Kcor_wm
#    print ggobs.weight, ggmw.weight, ggobs.npairs,ggmw.npairs

    # Write it out to a file and praise-be for Jarvis and his well documented code
    #gg.write(outfile)