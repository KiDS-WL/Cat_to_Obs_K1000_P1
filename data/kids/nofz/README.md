# SOM n(z) for KiDS-1000

Ascii histograms of the SOM-calibrated n(z) for each tomographic bin are located in the directory SOM_N_of_Z.    These are calculated following the methodology described in [Wright et al. 2020][1].

The correlated uncertainty on these n(z) is quantified through a covariance matrix for nofz shifts, calculated from 100 MICE2 LOS.  That covariance is then multiplied by 4 to take the uncertainity in the mocks into account. The results are in SOM_cov_multiplied.asc.   We also carry out a test analysis where we triple the error (SOM_cov_multiplied3.asc).

The MICE-calibrated offsets (deltaz) are given in deltaz.asc (see Table 3 in [Wright et al. 2020][1]). These are based on our best cuts with |<Z_B>-<z_spec>|<=0.4. 

In order to use correlated n(z) nuisance parameters we determine uncorrelated parameters using delta_z_correlated_to_uncorrelated.py.   This takes the covariance from SOM_cov_multiplied.asc and the delta_z from delta_z.asc and makes delta_z_unccorelated values through a cholesky decomposition of the covariance. We need delta_z_unccorelated values to be able to use the correlated_prior module in KCAP.

Finally we use a cross-correlation technique to validate the SOM n(z), as described in [Hildebrandt et al. 2020][2].  We use this analysis to update the priors on and correlation between the delta_z shifts.   The n(z) shape, however, remains fixed to the SOM n(z) values in the directory SOM_N_of_Z.

[1]: https://arxiv.org/pdf/1909.09632.pdf "Wright et al. 2020"
[2]: https://arxiv.org/abs/2007.15635 "Hildebrandt et al. 2020"
