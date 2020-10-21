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
The KiDS-1000 Cosmic Shear and 3x2pt covariance is described in [Joachimi, Lin, Asgari, Troester, Heymans et al. 2020][2].  The analytical covariance code will be made public shortly (please contact Benjamin Joachimi if you want early access). 
* covariance/inputs/blindC : initial ascii covariance matrix for the 2x2pt bandpower data vector, and the xi_pm cosmic shear data vector
* covariance/inputs/iterative_covariance/blindC : final ascii covariance matrix for the 2x2pt bandpower data vector, and the xi_pm cosmic shear data vector



[1]: https://arxiv.org/pdf/2007.15633.pdf "Asgari et al."





