import numpy as np

filename='deltaz.asc'
file=open(filename)
delta_z=np.loadtxt(file,comments='#')

filename='SOM_cov_multiplied.asc'
file=open(filename)
cov_z_in=np.loadtxt(file,comments='#')

cov_z = cov_z_in/4.*9.

filename='SOM_cov_multiplied3.asc'
np.savetxt(filename,cov_z,fmt='%.4e')


L = np.linalg.cholesky(cov_z) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z)

print(delta_x)


filename='SOM_cov_new.asc'
file=open(filename)
cov_z_in=np.loadtxt(file,comments='#')

cov_z = cov_z_in*4.


filename='SOM_cov_new_multiplied2.asc'
np.savetxt(filename,cov_z,fmt='%.4e')


L = np.linalg.cholesky(cov_z) 

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z)

print(delta_x)