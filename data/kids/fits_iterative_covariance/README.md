## Fits format KiDS-1000 data (final iterative covariance)

These fits cubes are in CosmoSIS format and include the tomographic data vector, covariance matrix and the n(z)'s.   The naming convention is as follows:

* bp = "band powers" - these are 2x2pt data vectors including the cosmic shear and GGL measurements
* cosebis = "COSEBIs" - these are cosmic shear only data vectors
* xipm = "two point correlation functions" - these are cosmic shear only data vectors
* bmode = "B-mode" - these contain the B-mode channel for band powers and COSEBIs
* with_m_bias = shear calibration correction applied to the data with the uncertainty carried through to the covariance matrix
* no_m_bias = shear calibration correction and uncertainty needs to be included as a nuisance parameter in the analysis

These files are created using 2pt_data_to_fits/save_and_check_Phase1.py

These files were used for the final published results using the iterated covariance matrix.
