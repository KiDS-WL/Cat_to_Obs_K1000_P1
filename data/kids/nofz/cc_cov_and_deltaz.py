import numpy as np


filename='deltaz_cc.asc'
file=open(filename)
delta_z=np.loadtxt(file,comments='#')


filename='CC_cov.asc'
file=open(filename)
cov=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z)

print('CC:',delta_x)



filename='CC_SOM_cov.asc'
file=open(filename)
cov=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z)

print('CC+SOM:',delta_x)


import scipy.linalg as la
import matplotlib.pylab as plt
filename='SOM_cov_multiplied.asc'
file=open(filename)
cov_som=np.loadtxt(file,comments='#')

# lambda_som, v = la.eig(cov_som)
# lambda_cc, v = la.eig(cov)

# plt.plot(np.sort(np.real(lambda_som)),'-r',marker='o')
# plt.plot(np.sort(np.real(lambda_cc)),'-b',marker='s')
# plt.show()
nBins=5
corr = cov.copy()
for i in range(nBins):
	for j in range(nBins):
		corr[i,j]= cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])

nBins=5
corr_som = cov_som.copy()
for i in range(nBins):
	for j in range(nBins):
		corr_som[i,j]= cov_som[i,j]/np.sqrt(cov_som[i,i]*cov_som[j,j])

vmin=min(corr.min(),corr_som.min()) 
vmax=max(corr.max(),corr_som.max())

fig, axes = plt.subplots(nrows=1, ncols=2)
im=axes[0].imshow(corr,vmin=vmin,vmax=vmax,aspect='auto')
axes[1].imshow(corr_som,vmin=vmin,vmax=vmax,aspect='auto')
plt.colorbar(im)
plt.show()

# plt.subplot(1,2,1)
# plt.imshow(corr)
# plt.subplot(1,2,2)
# plt.imshow(corr_som)
# plt.colorbar()
# plt.show()


