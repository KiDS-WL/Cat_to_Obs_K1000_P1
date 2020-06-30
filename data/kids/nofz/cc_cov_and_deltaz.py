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
