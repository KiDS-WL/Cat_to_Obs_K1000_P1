# Shear calibration corrections for the KiDS-1000 SOM sample

The KiDS-1000 shear estimates are calibrated through the analysis of image simulations that emulate KiDS, as described in [Kannawadi et al. 2019][1].  In KiDS-1000 we use the SOM n(z) as calibrated following the fiducial analysis of [Wright et al. 2020][2].    This imposes an object selection on the galaxies which we mimic in the image simulations to determine the calibration correction for the SOM sample (Summary_multiplicative_Fid_unblinded.npy, see Section 2.2 of [Giblin et al. 2020][3] for details).

For interest, we also provide the calibration corrections for other SOM samples, calibrated without DEEP2, VVDS or zCOSMOS, or requiring multiple different spectroscopic sources (see [Wright et al. 2020][2] for details).

We include this calibration correction by directly correcting the data vectors, and then including our uncertainty on the correction in the covariance matrix.  This can be done assuming correlated or uncorrelated but inflated errors.   This directory creates the relevant ascii files used to test how these choices impact our constraints (see Section 4.2.1 of [Asgari et al. 2020][4] for details).



[1]: https://arxiv.org/pdf/1812.03983.pdf "Kannawadi et al"
[2]: https://arxiv.org/abs/2005.04207 "Wright et al."
[3]: https://arxiv.org/abs/2007.01845 "Giblin et al."
[4]: https://arxiv.org/pdf/2007.15633.pdf "Asgari et al."
