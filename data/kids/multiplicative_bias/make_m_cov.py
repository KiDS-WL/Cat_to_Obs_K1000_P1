import numpy as np

# make m_cov based on sigma_m values for each redshift bin 
# and assume that there is 'corr' correlation between the bins. corr=1 is full correlation.
# m = np.asarray([-0.00930364, -0.01121353, -0.01518512,  0.00243974,  0.00754766])


sigma_m=np.asarray([0.019,0.020,0.017,0.012,0.010])
sigma_m_0p02=np.ones(len(sigma_m))*0.02

m_cov=np.diag(sigma_m**2)
corr = 0.99
str_corr = '%0.2f' % corr
nbins=len(sigma_m)
for i in range(nbins):
	for j in range(nbins):
		if i==j:
			print(m_cov[i,j])
		else:
			m_cov[i,j]=sigma_m[i]*sigma_m[j]*corr

np.savetxt('m_cov_r_'+str_corr+'.ascii',m_cov,fmt='%.4e')




m_cov=np.diag(sigma_m_0p02**2)
corr = 0.99
str_corr = '%0.2f' % corr
nbins=len(sigma_m)
for i in range(nbins):
	for j in range(nbins):
		if i==j:
			print(m_cov[i,j])
		else:
			m_cov[i,j]=sigma_m_0p02[i]*sigma_m_0p02[j]*corr

np.savetxt('m_cov_r_'+str_corr+'_0p02.ascii',m_cov,fmt='%.4e')


m_corr=np.zeros((nbins,nbins))
for i in range(nbins):
	for j in range(nbins):
		m_corr[i,j]=m_cov[i,j]/np.sqrt(m_cov[i,i]*m_cov[j,j])


# Here we make an uncorrelated m covariance which encompasses 
# all possible correlations that the original m_cov could have
# based on Hoyle et al.2018, appendix A: https://arxiv.org/pdf/1708.01532.pdf

m_cov = np.diag(len(sigma_m)*sigma_m**2)
np.savetxt('m_cov_uncorrelated_inflated.ascii',m_cov,fmt='%.4e')


m_cov = np.diag(len(sigma_m)*sigma_m_0p02**2)
np.savetxt('m_cov_uncorrelated_inflated_0p02.ascii',m_cov,fmt='%.4e')


