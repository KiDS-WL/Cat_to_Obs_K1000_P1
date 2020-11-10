# PSF Systematics:  quantifying and testing the impact of the Paulin-Henriksson and Bacon systematics model.
Codes written by B. Giblin

**Calc_PHterms.py**

This code calculates terms entering the Paulin-Henriksson systematics model given in Eqn eqn. 10 of Giblin et al. (2020).
It also calculates the error on the PH systematic contribution using Jackknife.

Arguments you can set in the code:
 - Calc_PH (True/False): whether to actually run TreeCorr & calculate the PH terms.
 - PreReadCat (True/False): whether to read in a pre-made pickled PSF catalogue to save time.
   If False it assembles the PSF catalogue from scratch using all the KiDS tile directories (slower).
   Always run it with PreReadCat set to	False the first	time, then True for subsequent runs.
 - Splitup_Fields (True/False): if True, splits observed field into Res*Res (RA,Dec) patches and uses
   these patches to calculate the Jackknife errors on the measurement.
 - Res (integer): For KiDS-1000, Res=7 (49 patches of ~20 deg^2 each seemed a good number).

 - LFver ("321"/"309b"): lensfit version used in producing the the PSF position-ellipticity catalogue.
   We used "321" for the PSF modelling
 - NorS ("N"/"S"): Use the data for the North or South KiDS fields.


**Plot_PHterms.py**

This code produces Figure 3 of Giblin et al. (2020) and also performs the chi^2 tests presented in that paper,
i.e. it checks that the induced change in chi^2 at a fiducial true cosmology that is caused by a given PSF
systematic to xi+, is some dominant to the change caused by a small increase or decrease in S_8.
Therefore it checks your PSF systematic does not bias your best-fit S_8 estimate overly.

Arguments:
 - Read_Tvalues (True/False): If False, this calculates the (delta)T_PSF/gal quantities appearing in, e.g.,
   the Paulin-Henriksson systematics model (eqn 10. Giblin+2020), from the data, per redshift bin.
   If True, it reads the pre-calculated values saved to file in order to save time on next run.
   Set to False on the first run and true on subsequent runs.
 - LFver: lensfit version. (see above).
 - Use_alpha_per_bin (True/False): alpha is the PSF leakage appearing in the Jarvis et al. (2016; eqn 3.16)
    systematics model. Note, we don't use this model in Giblin et al. (2020).
   If True, it reads an individual value of alpha for each redshift bin which were previously calculated using
   the Calc_1pt code.
   If False, it sets alpha to 0.03 for every bin (basically worst case scenario).

Notes:
  - To run the chi2 tests, execute the function Investigate_chi2.
  - To run these tests for different prescriptions for the PSF systematic (Paulin-Henriksson, Bacon, flux-dependent PSF residuals),
  - need to change the lines beginning "delta_xip... = " in this function (around line 337)
  - Calc_delta_xip_J16 = Jarvis+2016 systematic
  - Calc_delta_xip_H20 = Paulin-Henriksson+2008 systematic
  - Calc_delta_xip_cterms = flux-dependent PSF residual systematic
  - Calc_delta_xip_Bacon = Bacon+2003 systematic

**Marikas_CovMat/*.fits**

This directory contains two shear correlation function covariance matrices,
one including the multiplicative shear bias (*_with_m_bias_*)
and one without this component (*_no_m_bias_*).
This covariance is an analytical estimate described in detail in Appendix D of [Joachimi, Lin, Asgari, Tröster, Heymans et al. 2020][1].
This is the covariance used in the cosmic shear analysis of [Asgari, Lin, Joachimi et al. (2020)][2].

*STRUCTURE OF THE COVARIANCE:*

The covariance is measured in 9 logarithmically-spaced angular separation bins, saved under the 'ANG' column name.
The covariance matrix itself can be read in in python via, e.g.:
```python
from astropy.io import fits
filename = 'Marikas_CovMat/xipm_KIDS1000_BlindA_no_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits'
f = fits.open(filename)
cov = f[1].data
```

The covariance is shaped (270,270).
This shape corresponds to the 15 redshift bin combinations (5 auto-bins plus 10 cross-bin pairs),
times the 9 angular separation bins, times 2 as there are 2 shear correlation functions (xi+,xi-).
cov[0:135,0:135] is the xi+ covariance matrix, lowest redshift bin listed first.


**PSF_plots/average_delta_estar_x_y.py**

IF XY_or_chip is set to 'X_Y':
This code plots the mean & standard-deviation of the PSF ellipticity components (1&2), as well as the mean residuals, against the X,Y
position on the camera field-of-view.
i.e. it produces Figure 2 of Giblin et al. (2020).

IF XY_or_chip is set to 'chip':
it plots the mean PSF ellipticities and the residuals VS the r-band magnitude for each CCD separately,
i.e. it produces Figure 4 of Giblin et al. (2020).
The results are saved in PSF_plots/PLOTS/

This code uses PSF data that is produced by PSF_plots/paste_psf_residual_cats.sh
This code scrolls the PSF data catalogues saved for the various lines of sight creates a compiled data catalogue for use in the above code.



[1]: https://arxiv.org/abs/2007.01844 "Joachimi et al."
[2]: https://arxiv.org/abs/2007.15633 "Asgari et al."
