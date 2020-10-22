# Weighted galaxy pair counts

When calculating the shot-noise for the analytical covariance matrix we need to know the number of pairs in each bin.  This count needs to take into account the weights, however, as given in Appendix C.4 of [Joachimi, Lin, Asgari, Troester, Heymans et al. 2020][1].

In this directory you have the pair count for the broad 9-bin xi-files, compiled for all 15 tomographic bin combinations for cosmic shear e.g npair_blindC_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nTheta9.ascii,   and the fine binned GT and XI npairs, e.g pairCounts_K1000_blindC_Flag_SOM_Fid_NTheta326_nbTomo7.dat

[1]: https://arxiv.org/pdf/2007.01844.pdf "Joachimi et al."
