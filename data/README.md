# KiDS-1000 BOSS and 2dFLenS 3x2pt data

The KiDS-1000 cosmological analyses were blinded through the analysis of three different shear catalogues, Blind A, B and C.   Blind C was revealed as the truth after we completed our analysis.   In this repo to avoid any confusion we only provide the unblinded Blind C data, but you will find options in the scripts to create Blinds A and B.

## BOSS Sanchez et al 2017 Anisotropic Clustering
* boss/Sanchez_etal_2017:  clustering data and covariance matrix
* boss/nofz: BOSS and 2dFLenS n(z) for GGL predictions

## KiDS-1000 data
* kids/ellipticity_dispersion: sigma_e per tomographic bin
* kids/fits: fits data cubes for KCAP/CosmoSIS that contain the data vectors, initial covariance and n(z) as used in systematic tests
* kids/fits_iterative_covariance: fits data cubes for KCAP/CosmoSIS that contain the data vectors, final covariance and n(z) used for the final analysis
* kids/mock_data: mock data for pipeline tests
* kids/multiplicative_bias: m-calibration values from Kannawadi et al re-analysis of SOM sample
* kids/nofz: the SOM-calibrated n(z) and delta_z covariance matrix
* kids/number_density: neff for the KiDS-1000 area defined by the healpix nside 4096 mask used by the analytical covariance
* kids/number_of_galaxy_pairs: weighted pair counts for xi_pm (fine and broad) and gamma_t
* kids/psf_systematic_corrected: testing the impact of a PSF residual on the cosmological parameter inference
* kids/xipm: finely binned xi_\pm measurements per tomographic bin combination (ascii)
* kids also contains 15 tomographic-bin complication ascii files for the three statistics analysed in [Asgari, Lin, Joachimi et al. 2020][1], xi_\pm (broad), COSEBIs and band powers.

## Covariance Matrices
This directory holds the covariance matrices in ascii format.  If you want fits see the directory kids/fits.

The KiDS-1000 Cosmic Shear and 3x2pt covariance is described in Appendix E of [Joachimi, Lin, Asgari, Troester, Heymans et al. 2020][2].  This directory focuses on creating the COSEBIs covariance matrix, as described in Appendix A of [Asgari, Troester et al. 2020][3].  The analytical covariance code will be made public shortly (please contact Benjamin Joachimi and Marika Asgari if you want early access).   
* covariance/inputs/blindC : initial ascii covariance matrix for the 2x2pt bandpower data vector, and the xi_pm cosmic shear data vector
* covariance/inputs/iterative_covariance/blindC : final ascii covariance matrix for the 2x2pt bandpower data vector, and the xi_pm cosmic shear data vector
* covariance/inputs/input_for_xipm_sigma_m_covariance: theoretical xi_pm used to carry through the uncertainty on the shear calibration to the covariance matrix
* covariance/outputs: initial and final ascii covariance matrix for the COSEBIs cosmic shear data vector

## Data Plots
Maybe you don't want the data in CosmoSIS format?   No problem, take a look at the Data_Plots directory which contains the data and covariance in ascii format.   The python plotting scripts will guide you in how to incorporate the data into whatever format you want to work with it.  These scripts were used to create the figures in [Heymans, Tröster et al. 2020][4].

[1]: https://arxiv.org/pdf/2007.15633.pdf "Asgari et al. KiDS-1000"
[2]: https://arxiv.org/pdf/2007.01844.pdf "Joachimi et al."
[3]: https://arxiv.org/abs/1910.05336 "Asgari et al. KV450"
[4]: https://arxiv.org/abs/2007.15632 "Heymans et al."




