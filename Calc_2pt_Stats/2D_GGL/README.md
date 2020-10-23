# Code to measure the GGL of KiDS-1000 galaxies by BOSS lenses, in the OmegaCam pixel frame  (B. Giblin)

* Run_2D_GGL_K1000.py: 
Reads in KiDS-1000 shear catalogue and BOSS lens catalogue, applies designated redshift cuts to each and calculates
the 2D gamma_t signal in lots of patches out to 20 arcmin, stacks the patches.
Also can measure the 1D angular gamma_t with this same code.
This produces the measurements that go into Fig. 10 in Giblin et al. (2020). Using the 2D GGL signal to detect
additive shear systematics that vanish in the 1D angular GGL signal.

* Plot_2D_GGL_Results.py:
Plot the results from Run_GGL_K1000.py, producing Fig. 10 in Giblin et al. (2020).

* Functions_2DGGL.py:
Contains simple function that return elements in an array falling with a prescribed (ra1-ra2, dec1-dec2) patch.
Used by Run_2D_GGL_K1000.py.
