import numpy as np


# update this to the new one
filename="../thps_cov_kids1000_egretta_bp_apo_obs_bandpower_E_apod_list_mbias.dat"
file=open(filename)
m_bias_list=np.loadtxt(filename)[:,-1]

filename="thps_cov_kids1000_bandpower_E_apod_0_list.dat"
file=open(filename)
covariance_list=np.loadtxt(filename)

covariance_list_with_m = covariance_list.copy()
covariance_list_with_m[:,-1]+=m_bias_list

savename="thps_cov_kids1000_bandpower_E_apod_0_list_with_sigma_m.dat"
np.savetxt(savename,covariance_list_with_m,fmt='%i %i %i %i %i %i %i %i %i %i %1.9e')


