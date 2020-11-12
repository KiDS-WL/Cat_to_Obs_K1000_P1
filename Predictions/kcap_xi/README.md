This directory contains input files for kcap to produce xi+/- predictions at various values of S_8. These are read in by PSF_systests/Plot_PHterms.py to compare the impact of PSF-related biases on the xi+ to the changed caused by varying S_8.

 * **inputs/valueFiles/S8_<cosmol_tag>_values.ini**: input value files for kcap with S_8 varying slightly.
   - 'fid' is the fiducial S_8 value of 0.737.
   - All files have slight increases or decreases of S_8 by an amount determined by <cosmol_tag>. For example, high0.1sigma/low0.1sigma are increases/decreases of 0.1sigma_S8, where sigma_S8=0.020 is the approximate 1sigma error on S_8 in the KiDS-1000 cosmic shear analysis of Asgari et al. (2020).
   - The fiducial S_8 value, as well as the values of the other cosmological parameters in this ini file are the approximate best-fits of the KV450 cosmic shear analysis in Hildebrandt et al. (2018).
   - Outputs of this kcap script are saved in subdirectories labelled: Preditions/kcap_xi/outputs/test_output_S8_<cosmol_tag>_test/
   - In each subdirectory, there are further subdirectories for, e.g., the xi+/- predictions (/shear_xi_*), COSEBIs (/cosebis), band powers (bandpower_shear_e), the cosmological parameters (/cosmological_parameters), lin/nl matter power spectra (matter_power_*) etc.

 * **inputs/xi_variousS8.ini**: the input kcap file that specifies the pipeline. <cosmol_tag> can be set in this file to change which valuesFile is accessed.

 * **inputs/Run_all.ini **: runs the test sampler producing predictions for all two point statsitics. 
 
 * **inputs/nofz.fits **: a fits file that includes the redshift distribution of the galaxies. 

