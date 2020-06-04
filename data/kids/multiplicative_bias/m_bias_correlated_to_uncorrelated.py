import numpy as np


# Here we find the values for the centre of the Gaussian priors for the uncorrelated m_bias parameters

filename="Summary_multiplicative_Fid_unblinded.npy"
m=np.load(filename)[:,1]

filename='m_cov_r_0.99.ascii'
file=open(filename)
cov_m=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov_m) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,m)

print('the correlated case:', delta_x)


filename='m_cov_r_0.99_0p02.ascii'
file=open(filename)
cov_m=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov_m) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,m)

print('the correlated case with 0.02 sigma_m:', delta_x)



filename='m_cov_uncorrelated_inflated.ascii'
file=open(filename)
cov_m=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov_m) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,m)

print('the uncorrelated case:',delta_x)



filename='m_cov_uncorrelated_inflated_0p02.ascii'
file=open(filename)
cov_m=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov_m) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,m)

print('the uncorrelated case with 0.02 sigma_m:',delta_x)



