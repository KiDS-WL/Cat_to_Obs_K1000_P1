21/10/2020: codes written by H. Hildebrandt & B. Giblin, hacked by B.Giblin and C.Heymans.

1. create_tomocats_GGL.py:
For the shear ratio test we need finer lens bins than normal so we have a special script to create these

2. run_shear_ratio_w_spin.sh:
   For a given source bin (specified on command line), this code runs Shear_ratio_wspin_test.py recursively for several lens bins.
   Accepts a parameter file which specifies whether to run this for KiDS-1000 or MICE mocks (more on that below).


3. Shear_ratio_wspin_test.py:  the workhorse code that makes all the measurements
Calculate the gamma_t measurement between source and lens samples specified in the input parameter file.
It also	calculates the predictions for the covariance of the gamma_t via different methods
(in Giblin+2020 we use spin-shear realisations).	
The TreeCorr memory leak should	now be fixed.
 Arguments:
 - Cov_Method: option to change	the covariance being calculated	from multiple spun-shear realistions
   to instead read in the octant-sized MICE simulation and split this into n patches (with n~16).
 - DFLAG: only used if running test with MICE mocks. If	set to '_Octant' it will use the MICE octant-sized simulation
   as the data catalogue in the	test.
   Default is to leave DFLAG='' which uses the 343 sqdeg KV450-like MICE catalogue.
 (In input parameter file):
 - SOURCE_TYPE:	change to using	the KiDS-1000 or MICE mock catalogues for the sources.
 - LENS_TYPE: same as above but	for the	lens sample (with MICE lenses, can have BOSS-like or GAMA-like lenses).
 - Mag_OnOff: turn the magnification on/off in the MICE mocks.
 - Pz_TrueEstimated: whether the source redshift n(z) is the truth or the estimated in the MICE mocks.
 - SN: galaxy shape noise on or off.


4. calc_spin_test_cov.py - creates a covariance matrix from all the spin trials
Once the predictions for the covariance estimation have been made by run_shear_ratio_w_spin.sh / Shear_ratio_wspin_test.py,
this code reads them in, calculates the covariane matrix and saves it.

Once the data is assembled the main shear ratio test is done in

5. GGL_shear_ratio_test_zall.py  -
This makes Figure 11 of Giblin et al. (2020).

this needs some input data however from...

6. Dls_over_Ds.py: which calculates a bunch of different shear ratios for a fixed cosmology (using cosmology.py)

7. run_Dls_over_Ds.sh: is the master script that creates the Dls/Ds files

8. Compare_BFParams_And_Models.py
   This code was used to read in the best-fit gamma_t model and corresponding B_ij model parameters (eqn 18 Giblin+2020)
   and compare them for different cases. These cases included:
   - n(z) systematically shifted up or down (or incoherently shifted, up/down depending on zbin).
   - The params/models fitted when each redshift bin was modelled separately.
   - having a z=1.4 peak artifically added to the n(z) to see the effect of potential high redshift outliers that are unconstrained
     by our photo-z calibration.
   - performing the SRT for all 5 source bins, or just the 3 highest source bins.
   - Including/neglecting the uncertainty due to the intrinsic alignment in the covariance.
   - whether an extra parameter per lens bin was included to model the contribution of magnification to the gamma_t
   The results of these tests are discussed to some extent in Sect. 4.3.2 of Giblin+2020.
   Concluded that the SRT is highly insensitive to these calibration errors given the statistical power of Stage-III surveys,
   especially when IA-uncertainty is included.
   Conclusion: SRT is not a very good test.

9. Make_MockMagnifiedgt.py
  This was used to artificially add a magnification contribution to a mock gamma_t, with a free parameter per lens bin,
  alpha, dictating the strength of the added contribution.
  This was used to check that the code, GGL_shear_ratio_test_zall.py, could correctly fit the alpha-per-lens bin
  when the option to include extra parameters to model magnification was turned on in that code. 