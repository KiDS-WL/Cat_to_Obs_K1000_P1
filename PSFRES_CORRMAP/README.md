# Chip-dependent 2D PSF residual model

In Section 3.4.1 of [Giblin, Heymans, Asgari et al. 2020][1] we discuss the flux-dependent detector-level distortions which lead to a chip-dependent PSF residual error.   The baseline PSF residual model is shown in Figure 2 of [Hildebrandt et al. 2020][2], with the negligible impact on the expected cosmological signal shown in Fig 7 of [Giblin, Heymans, Asgari et al. 2020][1].

This directory contains the PSF residual model extrapolated to faint magnitudes and dithered to form a co-add image (c1_map.fits and c2_map.fits) and code to create a mock KiDS catalogue of this PSF residual model at the location of each KiDS galaxy (create_c12_mock.py).   Scripts are also provided to fit a free amplitude for the model directly to the data.  

[1]: https://arxiv.org/abs/2007.01845 "Giblin et al."
[2]: https://arxiv.org/abs/1812.06076 "Hildebrandt et al."
