This directory contains input files for kcap to produce xi+/- predictions at various values of S_8. These are read in by PSF_systests/Plot_PHterms.py to compare the impact of PSF-related biases on the xi+ to the changed caused by varying S_8.

 * **inputs/valueFiles/S8_<cosmol_tag>_values.ini**: input value files for kcap with S_8 varying slightly. 'fid' is the fiducial S_8 values. All others are slight increases or decreases by an amount specified by <cosmol_tag>.

 * **inputs/xi_variousS8.ini**: the input kcap file that specifies the pipeline. <cosmol_tag> can be set in this file to change which valuesFile is accessed.

 * **inputs/Run_all.ini **: not sure what this is - I think it is from Marika...? It uses the nofz.fits available in this directory.

