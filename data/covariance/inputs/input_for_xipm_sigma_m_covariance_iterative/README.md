# CosmoSIS config files and outputs used for getting theoretical xipm to be used in the covariance calculations. 

The uncertainity in the multiplicative bias is added to the covariance matrix in the fiducial analyses of KiDS-1000, 
following equation 37 of [Joachimi, Lin, Asgari, Troester, Heymans et al. (2020)][1].
We assume that this uncertainity is fully correlated between the tomographic bins and use a Gaussian with a width estimated from [Kannawadi et al. (2019)][2]. 
See the last column of table 1 in [Asgari, Lin, Joachimi et al. (2020)][3] for the values of sigma_m. 

The python script: add_sigma_m_cov_to_xipm.py in ../blindC/ takes the theoretical xipm saved in chains/ and adds the sigma_m contribution to the xipm covariance matrix.

The input_for_xipm_sigma_m_covariance directory contains the files for our initial analysis using a fiducial cosmology taken from the KV450+BOSS analysis of [Troester et al. 2020][4].   This directory updates the input cosmology to the best-fit cosmology from the first 3x2pt analysis with the fiducial cosmology.


[1]: https://arxiv.org/abs/2007.01844 "Joachimi et al."
[2]: https://arxiv.org/abs/1812.03983 "Kannawadi et al."
[3]: https://arxiv.org/abs/2007.15633 "Asgari et al."
[4]: https://arxiv.org/abs/1909.11006 "Troester et al."
