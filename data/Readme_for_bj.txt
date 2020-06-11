The data you need is in the kids folder.

sigma_e in kids/ellipticity_dispersion
n(z)    in kids/nofz/SOM_N_of_Z
n_eff   in kids/number_density
nPairs  in kids/number_of_galaxy_pairs

For nPairs you have the pair count for the broad 9-bin xi-files e.g npair_blindA_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nTheta9.ascii,   and the fine binned GT and XI npairs, e.g pairCounts_K1000_blindA_Flag_SOM_Fid_NTheta326_nbTomo7.dat



If you want the sigma_e and neff in the format that you had for Blind A, you can find K1000_Phase1_ZB_bin_sigma_e_neff_blind*_FLAG_SOM_DR4.1.txt in the kids/ellipticity_dispersion folder

Please put the results in 

data/covariance/inputs/blindB 
and 
data/covariance/inputs/blindC

There is already a blindA folder with the inputs in there. 
