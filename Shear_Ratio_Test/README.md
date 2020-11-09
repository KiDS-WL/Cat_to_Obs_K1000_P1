# Shear Ratio Test

The code in this directory performs the KiDS-1000 shear ratio test described in Section 4.3 of [Giblin, Heymans, Asgari et al. 2020][1].   The code was originally written by H. Hildebrandt & B. Giblin, and hacked by B. Giblin and C. Heymans.

## Creating Lens and Source catalogues

**@Ben - Add code/directory description here - e.g what are the LENSCAT/SOURCECAT directories?**  

* create_tomocats_GGL.py:
For the shear ratio test we need finer lens bins than normal so we have a special script to create these


## The Shear Ratio Test

* **Shear_ratio_wspin_test.py**:  the workhorse code that makes all the measurements
Calculate the gamma_t measurement between source and lens samples specified in the input parameter file.
It also	calculates the predictions for the covariance of the gamma_t via different methods
(in Giblin+2020 we use spin-shear realisations).	

 Arguments:
   - Cov_Method: option to change	the covariance being calculated	from multiple spun-shear realistions
   to instead read in the octant-sized MICE simulation and split this into n patches (with n~16).
   - DFLAG: only used if running test with MICE mocks. If	set to '_Octant' it will use the MICE octant-sized simulation
   as the data catalogue in the	test.
   Default is to leave DFLAG='' which uses the 343 sqdeg KV450-like MICE catalogue.
   - (In input parameter file):
    - SOURCE_TYPE:	change to using the KiDS-1000 or MICE mock catalogues for the sources.
    - LENS_TYPE: same as above but for the	lens sample (with MICE lenses, can have BOSS-like or GAMA-like lenses).
    - Mag_OnOff: turn the magnification on/off in the MICE mocks.
    - Pz_TrueEstimated: whether the source redshift n(z) is the truth or the estimated in the MICE mocks.
    - SN: galaxy shape noise on or off.
    
* **run_shear_ratio_w_spin.sh**:
   For a given source bin (specified on command line), this code runs Shear_ratio_wspin_test.py recursively for several lens bins.
   Accepts a parameter file which specifies whether to run this for KiDS-1000 or MICE mocks (more on that below).

## Calculating the Covariance

* **calc_spin_test_cov.py**: creates a covariance matrix from all the spin trials.
Once the predictions for the covariance estimation have been made by run_shear_ratio_w_spin.sh / Shear_ratio_wspin_test.py,
this code reads them in, calculates the covariance matrix and saves it.

## Calculating the expected redshift-scaling

* **Dls_over_Ds.py**: given a lens and source n(z), this code calculates the shear ratio for a fixed cosmology (using cosmology.py).   Options are included to shift the source n(z) by an error +/- delta_z, and also to include and n(z) with an artifical high-z outlier population included.
* **run_Dls_over_Ds.sh**: is the master script that creates the Dls/Ds files, running over all variant combinations of the source and lens bins.

## ADD SUITABLE TITE
* **GGL_shear_ratio_test_zall.py**  -
This makes Figure 11 of Giblin et al. (2020).
It does this by reading in the gamma_t^ij measurements corresponding to the i'th/j'th lens/source bin pairs,
reading in the covariance produced by calc_spin_test_cov.py,
reading in the Dls/Ds luminosity distance ratios (ls=lens-to-source, s=source) prouced by Dls_over_Ds.py,
and it fits a model given by eqns 18 and 19 of Giblin et al. (2020).
The model has one free amplitude parameter per lens bin.
On top of this, there's loads of extra options you can toggle. These include:
- Include_mCov (True/False)            # Include the uncertainty due to the m-correction
- Include_mBias (True/False)           # If True, read in and bias the gamma_t by *(1 + m + f_mBias * m_err)
                                       # where f_mBias, the level of bias to apply, can be set in the code).
- Include_Hartlap (True/False)         # Inflate the covariance by the Hartlap scaling factor, to account for noise in taking inverse.
- Include_Magnification (True/False)   # If True, include extra param in gamma_t model: strength of magnifcation effect on gt.
- Include_IA = True                    # If True, read in the kcap prediction for the IA-only <g_t> and inflate the cov diag by
- f_IA = 0                             # +(f_IA* gt_IA)^2 ; f_IA is our uncert. on IA amplitude.
- A_IA = 1                             # A_IA*gt_IA is added to the fitted model if Include_IA is True.

- nofz_shift="_ModnofzUpDown"          # Only for K1000: use the Dls/Ds values for the nofz which has been
                                       # shifted up ('_nofzUp'), down ('_nofzDown') by +/-(delta_z+delta_z_err)
                                       # OR 5sig shift-up ('_nofzUp5sig'), down ('_nofzDown5sig') by (+/- 5*delta_z_err)
                                       # "_nofzMix5sig" means alternate bins are shifted up/down by (+/- 5*delta_z_err).
                                       # For no shift, set to ''.
                                       # Finally, to include the uncert. on the nofz's in the SRT modelling,
                                       # set nofz_shift to "_ModnofzUpDown"
- Cov_Method="Spin"/"Patch"            # Spin means the cov was calculated from many gamma_t realisations
                                       # with the source galaxy ellipticities randomised.
                                       # "Patch" means the cov came from the MICE octant-sized simulation, split into ~16 patches.

GGL_shear_ratio_test_zall.py gets its input data from the following codes:


## Magnification and Intrinsic Galaxy Alignments

* **Compare_BFParams_And_Models.py**:
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

* **Make_MockMagnifiedgt.py**:
  This was used to artificially add a magnification contribution to a mock gamma_t, with a free parameter per lens bin,
  alpha, dictating the strength of the added contribution.
  This was used to check that the code, GGL_shear_ratio_test_zall.py, could correctly fit the alpha-per-lens bin
  when the option to include extra parameters to model magnification was turned on in that code. 
  
  
  
  [1]: https://arxiv.org/pdf/2007.01845.pdf "Giblin et al."
