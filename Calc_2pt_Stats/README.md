# Calculating 2pt Statistics from the KiDS-1000 shear and photo-z catalogues
This directory is the main workhouse of the repo, taking the KiDS-1000 shear and photo-z catalogues and converting them into the 2pt statistics that are then used in the cosmological parameter inference.   The main workhorse for the scripts is the fantastic [TreeCorr][1] code from Mike Jarvis.

## doall_calc2pt.sh
We're from the school of Thomas Erben who likes doall scripts!   This script really "does all".   If you like bash scripts, you'll be happy.   If you don't like bash scripts, don't worry, you don't need to use this, but you can read it to see how the different python/c-code can be called with the relevant command line information.   The doall script has the following modes
* CREATETOMO : uses create_tomocats.py to cut catalogues into tomographic bins and calculate and subtract c-term
* GT: uses calc_gt_w_treecorr.py to calculate gamma_t/x for a given tomographic bin combination in the North or South patch
* XI: uses calc_xi_w_treecorr.py to calculate xi_pm for a given tomographic bin combination in the North or South patch
* COMBINEGT: combines the gamma_t/x results from KiDS North/South for a given tomographic bin combination
* COMBINEXI: combines the xi_pm results from KiDS North/South for a given tomographic bin combination
* REBINGT: we calculate gamma_t/x in fine bins for the bandpower calculation, this mode uses calc_rebin_gt_xi.py to rebin the fine data into broad bins
* REBINXI: we calculate xi_pm in fine bins for the bandpower calculation, this mode uses calc_rebin_gt_xi.py to rebin the fine data into broad bins
* Pgk: uses src/bandpowers/xi2bandpow.c to convert gamma_t into P-galaxy-matter bandpowers
* Pkk: uses src/bandpowers/xi2bandpow.c to convert xi_pm into P-matter-matter bandpowers
* COSEBIs: uses src/cosebis/run_measure_cosebis_cats2stats.py to convert xi_pm into COSEBIs modes

The doall_calc2pt.sh script is called by cuillinrun_calc2pt.sh which, for the most intensive tasks, uses slurm to send jobs to different workers on cuillin (the ROE compute cluster).  We also have cuillinrun_calc2pt_mocks.sh which is the queue scheduler for the SALMO mocks.

The list of modes above cover most of the scripts in this directory, but we have some additions including
* create_tomocats_GGL.py, which creates finely z-binned BOSS catalogues for the shear ratio test
* create_tomocats_with_depsf_XY.py with the delta_epsf(X,Y) model from the directory PSFRES_CORRMAP subtracted
* ldac.py:  KiDS catalogues are fits tables and ldac-compatible.   You can use astropy functions to access the fits tables, but in most instances ldac.py is faster.

## 2D_GGL
In Section 4.2 of [Giblin, Heymans, Asgari et al. 2020][2] we present galaxy-galaxy lensing measurements in the OmegaCAM pixel reference frame.   Run_2D_GGL_K1000.py uses [TreeCorr][1] with [bin_type = "TwoD"][3] to carry out this measurement.  Plotting scripts are also provided in this directory.

## PLOTTING
This directory contains simple plotting scripts used to create Figures in [Giblin, Heymans, Asgari et al. 2020][2].
* plot_COSEBIs_Bmodes.py:  Figure 9
* plot_CSys.py: Figure 8
* plot_bandpower_Bmodes.py: Figure not included in Giblin et al to save space
* plot_star_gal_comp.py: Figure 7



[1]: https://github.com/rmjarvis/TreeCorr "TreeCorr from Mike Jarvis"
[2]: https://arxiv.org/pdf/2007.01845.pdf "Giblin et al."
[3]: https://rmjarvis.github.io/TreeCorr/_build/html/binning.html#twod "bin_type = TwoD"
