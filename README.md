# Cat_to_Obs_K1000_P1

This repository contains the scripts used to convert the shear-photo-z catalogues from the fourth data release of the Kilo-Degree Survey (KiDS-1000) into 1pt and 2pt observables.  The software was used in the following publications

* [Giblin, Heymans, Asgari et al. 2020][1]: _KiDS Catalogues: weak gravitational lensing shear measurements_
* [Joachimi, Lin, Asgari, Tröster, Heymans et al. 2020][2]: _KiDS-1000 Methodology: Modelling and inference for joint weak gravitational lensing and spectroscopic galaxy clustering analysis_
* [Asgari, Lin, Joachimi et al. 2020][3]: _KiDS-1000 Cosmology: Cosmic shear constraints and comparison between two point statistics_
* [Heymans, Tröster et al. 2020][4]: _KiDS-1000 Cosmology: Multi-probe weak gravitational lensing and spectroscopic galaxy clustering constraints_

With the publication of these papers, we make these scripts available to the community following the draft Software Policy document from the LSST-DESC Software Review Policy Committee.   In doing so we aim to show the reliability, reproducibility and reusability of our code, ensure visibility for our developers, and build trust in our "4-eye" software review process.   There is a major caveat here though - we are not software engineers! We fully recognise that our code could be improved in many ways.  We nevertheless provide these scripts "as is" without user support.

We provide a brief summary of each directory, referring the user to the directory readmes which further describe the contents.  
* 2pt_data_to_fits: Creating data fits tables for KCAP/CosmoSIS
* Calc_1pt_Stats: Calculating 1pt Statistics from the KiDS-1000 shear and photo-z catalogues
* Calc_2pt_Stats: Calculating 2pt Statistics from the KiDS-1000 shear and photo-z catalogues
* Flinc_theory_inputs: [SALMO][5] mock inputs
* GGL_LensCats: Creating KiDS-overlap Lens Catalogues from the BOSS and 2dFLenS master catalogues
* PSFRES_CORRMAP : Chip-dependent 2D PSF residual model
* PSF_systests: PSF Systematics: quantifying and testing the impact of the Paulin-Henriksson and Bacon systematics model.
* data: KiDS-1000 BOSS and 2dFLenS 3x2pt data and covariances in ascii and CosmoSIS-friendly fits format
* src: Code to convert finely binned xi_\pm measurements into COSEBIs and Bandpower estimates







[1]: https://arxiv.org/abs/2007.01845 "Giblin et al."
[2]: https://arxiv.org/abs/2007.01844 "Joachimi et al."
[3]: https://arxiv.org/abs/2007.15633 "Asgari et al."
[4]: https://arxiv.org/abs/2007.15632 "Heymans et al."
[5]: https://github.com/Linc-tw/salmo "SALMO"
