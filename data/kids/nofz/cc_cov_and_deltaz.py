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


# import scipy.linalg as la
# import matplotlib.pylab as plt
# filename='SOM_cov_multiplied.asc'
# file=open(filename)
# cov_som=np.loadtxt(file,comments='#')

# lambda_som, v = la.eig(cov_som)
# lambda_cc, v = la.eig(cov)

# plt.plot(np.sort(np.real(lambda_som)),'-r',marker='o')
# plt.plot(np.sort(np.real(lambda_cc)),'-b',marker='s')
# plt.show()

# plt.subplot(1,2,1)
# plt.imshow(cov_som)
# plt.subplot(1,2,2)
# plt.imshow(cov)
# plt.show()


