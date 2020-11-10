import numpy as np

Foldername = 'Cat_to_Obs_K1000_P1/data/kids/nofz/'

# specify the bin(s) to remove in a list
bins_to_remove=[5] # options are [1], [2], [3], [4], [5], [1,2]

str_bin_names=''
for bin_i in bins_to_remove:
	str_bin_names+=str(bin_i)
	str_bin_names+='_'

# Read in the som covariance --> use this coavriance instead of the fiducial ones in the pipeline.ini file
file=open(Foldername+'SOM_cov_multiplied_bin'+str_bin_names+'removed.asc')
cov=np.loadtxt(file)

# get the cholesky decomposition
L = np.linalg.cholesky(cov) 
inv_L = np.linalg.inv(L)

# read in the delta_z values for the SOM dist caluclated from MICE
filename=Foldername+'deltaz.asc'
file=open(filename)
delta_z_in=np.loadtxt(file,comments='#')

# Use these values for the uncorrelated parameters. 
delta_x = np.dot(inv_L,delta_z_in)

delta_z_out = np.dot(L,delta_x)

print('bins_to_remove:',bins_to_remove)
print('delta_z_in=',delta_z_in)
print('delta_x=',delta_x)
print('delta_z_out=',delta_z_out)


