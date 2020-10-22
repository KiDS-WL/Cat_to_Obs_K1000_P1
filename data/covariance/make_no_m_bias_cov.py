import numpy as np
import matplotlib.pylab as pl



blind = 'C'
filename='outputs/Covariance_bestfit_3x2pt_blind'+blind+'_nMaximum_20_0.50_300.00_nBins5.ascii'
file=open(filename)
cov_all=np.loadtxt(file,comments='#')

filename='outputs/Covariance_bestfit_3x2pt_blind'+blind+'_nMaximum_20_0.50_300.00_nBins5_sigma_m_from_m_cov_0.0000.ascii'
file=open(filename)
cov_sigma_m=np.loadtxt(file,comments='#')

cov_no_m_bias = cov_all - cov_sigma_m

filename='outputs/Covariance_bestfit_3x2pt_no_m_bias_blind'+blind+'_nMaximum_20_0.50_300.00_nBins5.ascii'
np.savetxt(filename,cov_no_m_bias,fmt='%.20e')


# filename='outputs/Covariance_no_m_bias_blindA_nMaximum_20_0.50_300.00_nBins5.ascii'
# file=open(filename)
# cov_no_m_bias_a = np.loadtxt(file,comments='#')

# pl.imshow((cov_no_m_bias/cov_no_m_bias_a)-1)
# pl.colorbar()
# pl.show()
