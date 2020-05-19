import numpy as np

filename='deltaz.asc'
file=open(filename)
delta_z=np.loadtxt(file,comments='#')



filename='SOM_cov_multiplied2.asc'

file=open(filename)
cov_z=np.loadtxt(file,comments='#')

L = np.linalg.cholesky(cov_z) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z)

print(delta_x)