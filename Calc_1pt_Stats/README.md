# Calculating 1pt Statistics from the KiDS-1000 shear and photo-z catalogues

Codes written primarily by B. Giblin, hacked by C. Heymans.

1. Plot_1pt_Tests.py:
This code reads in the shear catalogue and the PSF data catalogue for KiDS-1000 and produces plots including:
   - avg PSF ellipticity and residual PSF ellipticity VS (X,Y) position on CCD (Fig. 2, Giblin et al. 2020).
   - avg galaxy ellipticity VS (X,Y). 
   - avg y VS avg x : e.g., alpha (PSF leakage) VS tomo redshift bin (Fig. 6 of Giblin+20). 
   - avg magnitude VS ZB
 etc.

2. Plot_Star_Col-Col.py:
This code reads in the KiDS-1000 source catalogue & PSF catalogue and plots the colour-colour distribution of these objects
(Fig. 1 of Giblin+20). By cross-matching the objects in the source and PSF catalogues, we quantify the galaxy contamination
to the PSF sample. 
Overplots the Baldri et al. (2010) stellar locus.

3. fitting.py:
Simple code to fit either a linear function, or a general user-defined function, to some input data.

4. Make_K1000_MassMaps.py:
Read in shear data for KiDS-1000 and perform spherical harmonic transform to get kappa maps.
Also do a series of noise realisations to make signal/noise maps.
The products of this code are not intended for scientific analysis, mainly just for pretty maps / press release.

5. K1000_MassMaps_SimplePlotCode_4KiDS.py:
Read in and plot different projections of the KiDS-1000 SNR maps (kappa maps / noise maps).
This code was made to optimise the prettiness of the maps.
