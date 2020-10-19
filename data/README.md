# KiDS-1000 BOSS and 2dFLenS 3x2pt data

The KiDS-1000 cosmological analyses were blinded through the analysis of three different shear catalogues, Blind A, B and C.   Blind C was revealed as the truth after we completed our analysis.   In this repo to avoid any confusion we only provide the unblinded Blind C data, but you will find options in the scripts to create Blinds A and B.

## BOSS Sanchez et al 2017 Anisotropic Clustering
* boss/Sanchez_etal_2017:  clustering data and covariance matrix
* boss/nofz: BOSS and 2dFLenS n(z) for GGL predictions

## KiDS-1000 data
* kids/cosebis: Ascii E and B mode COSEBIS for each tomographic bin combination (*to do add Blind C*)
* kids/ellipticity_dispersion: sigma_e per tomographic bin
* kids/fits: fits data cubes for KCAP/CosmoSIS that contain the data vectors, initial covariance and n(z) as used in systematic tests
* kids/fits_iterative_covariance: fits data cubes for KCAP/CosmoSIS that contain the data vectors, final covariance and n(z) used for the final analysis
* kids/mock_data: mock data for pipeline tests
* kids/multiplicative_bias: m-calibration values from Kannawadi et al re-analysis of SOM sample



* sigma_e in kids/ellipticity_dispersion
* n(z)    in kids/nofz/SOM_N_of_Z
* n_eff   in kids/number_density
* nPairs  in kids/number_of_galaxy_pairs

For nPairs you have the pair count for the broad 9-bin xi-files e.g npair_blindA_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nTheta9.ascii,   and the fine binned GT and XI npairs, e.g pairCounts_K1000_blindA_Flag_SOM_Fid_NTheta326_nbTomo7.dat



