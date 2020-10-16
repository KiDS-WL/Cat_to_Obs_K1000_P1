# Creating KiDS-overlap Lens Catalogues from the BOSS and 2dFLenS master catalogues

Author:  Chris Blake and Catherine Heymans

Questions to: cblake@swin.edu.au

The code uses: numpy, scipy, astropy, matplotlib.

## Download BOSS/2dFLenS master catalogue data

The original datasets you need:

(1) The BOSS DR12 LRG data and random catalogue in the NGP:

* https://data.sdss.org/sas/dr12/boss/lss/
* wget https://data.sdss.org/sas/dr12/boss/lss/galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz
* wget https://data.sdss.org/sas/dr12/boss/lss/random0_DR12v5_CMASSLOWZTOT_North.fits.gz
* wget https://data.sdss.org/sas/dr12/boss/lss/random1_DR12v5_CMASSLOWZTOT_North.fits.gz

(2) The 2dFLenS LRG data and random catalogue in both KIDS-S and KiDS-N, I have prepared new versions on the 
2dFLenS website in the revised DR12 redshift slices here (http://2dflens.swin.edu.au).  Click on "data release" -- the 6 catalogues under "data generated in final BOSS redshift binning"

* wget http://2dflens.swin.edu.au/data_2dfbz1_kidss_160105_rat.tar.gz
* wget http://2dflens.swin.edu.au/data_2dfbz1_kidsn_160105_rat.tar.gz
* wget http://2dflens.swin.edu.au/data_2dfbz2_kidss_160105_rat.tar.gz
* wget http://2dflens.swin.edu.au/data_2dfbz2_kidsn_160105_rat.tar.gz
* wget http://2dflens.swin.edu.au/data_2dfbz3_kidss_160105_rat.tar.gz
* wget http://2dflens.swin.edu.au/data_2dfbz3_kidsn_160105_rat.tar.gz

(3) The KiDS-S and N superuser photometry catalogues. Note these need to be superuser KiDS catalogues, not the lensing catalogues, as the lensing catalogues have a weight>0 cut which removes most of the BOSS/2dFLenS galaxies

## makelenscats.py

The makelenscats.py code performs the following steps:

(1) Specify redshift bin to generate (ired = 1/2 on the command line).

(2) Read in BOSS data/randoms in KiDS-N, cutting to the KiDS-1000 footprint as specified by the KiDS-BOSS-2dFLenS Healpix Mask

(3) Read in 2dFLenS data/randoms in KiDS-S 

(4) Read in KiDS photometric catalogues in KiDS-S and KiDS-N.

(5) Match BOSS/2dFLenS data to KiDS photometric catalogues, nearest neighbour within 1 arcsec.  
    There is a flag "iphot" which tracks whether or not there is a match.  
    The photometric quantities extracted for the magnitude/colour weighting are: g-r, r-i and r. 

(6) Determine weights to match magnitude/colour distribution of 2dFLenS to BOSS, 
    using Hendrik's KV450 DIR method.  Unmatched objects are given weight=1.

(7) Determine weights to match magnitude/colour distribution of BOSS to 2dFLenS.  
    Unmatched objects are given weight=1.

The code then outputs fits files for all the data and random lenses (not just matched lenses) with columns:

RA,DEC,Z -- R.A., Dec., redshift
WEICOMP -- BOSS completeness weight (=1 for 2dFLenS)
WEIFKP -- FKP weights (these are used for BOSS clustering, so you probably need them)
FLAGPHOT -- 0 or 1 if there is a match to the photometric catalogues
WEIMAG -- weight to match the magnitude/colour distribution of the other survey
GRCOL, RICOL, RMAG -- saving these in order to easily produce plots of the distributions 
                      (could potentially delete these columns)

Random lenses have WEICOMP=1, FLAGPHOT=0, WEIMAG=1 
                      (not entirely sure what weight to give random lenses if the data lenses have 
                      magnitude/colour weights).

The standard weight to use to correspond to the BOSS clustering analysis would be WEICOMP*WEIFKP.

The other branch of the code ("testcats") produces various check plots and outputs which
produces a useful weighted vs unweighted plot of the colours and magnitudes amongst other things

## Other code

* BOSS_2dFLenS_n_of_z_calc_and_plot.py does what it says on the tin!
* The CREATE_FITS_HEALPIX_MASKS directory contains code to create the mosaic healpix mask for 2dFLenS and BOSS.   There are also simple scripts to calculate the overlap area and the effective number density of lenses in each redshift bin.
