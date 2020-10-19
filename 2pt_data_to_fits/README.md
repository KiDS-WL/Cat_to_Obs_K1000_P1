# Creating fits tables for KCAP
2pt data (xi_pm, COSEBIS, bandpower, gamma_t) from is created using software in the Calc_2pt_Stats directory, and stored as ascii output files for individual tomographic bin combinations.   This directory contains python scripts to convert these ascii outputs to a single fits table, that also includes the n(z) and covariance, for analysis with [KCAP][1].   

## save_and_check_Phase1*.py
The main workhorse that takes the relevant input data files, applies any scale cuts, and outputs the fits table in CosmoSIS-friendly format. The "iteration.py" version was used to make the final KiDS-1000 fits table with the iterative covariance matrix calculated at the best-fit cosmology of the initial run.  The initial run was used for all the systematic tests. We therefore provide both versions, where the main difference is the name of some of the input files.

## MakeDataVectors*.py
If you only want to combine the different tomographic bins use this code.   It also allows for the rebinning of the finely measured xi_pm results into broad bins

[1]: https://github.com/KiDS-WL/kcap "KCAP"
