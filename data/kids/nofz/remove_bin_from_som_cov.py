import numpy as np

file=open('SOM_cov_multiplied.asc')
cov_in=np.loadtxt(file)

L_in = np.linalg.cholesky(cov_in) 
inv_L_in = np.linalg.inv(L_in)

bins_to_remove=[5]

# remove row and column from cov
cov = cov_in.copy()
str_bin_names=''
for bin_i in bins_to_remove:
	cov[bin_i-1,:]=0.
	cov[:,bin_i-1]=0.
	cov[bin_i-1,bin_i-1]=1.
	str_bin_names+=str(bin_i)
	str_bin_names+='_'


filename='SOM_cov_multiplied_bin'+str_bin_names+'removed.asc'
np.savetxt(filename,cov,fmt='%.4e')


L = np.linalg.cholesky(cov) 

filename='deltaz.asc'
file=open(filename)
delta_z_in=np.loadtxt(file,comments='#')

delta_x_in = np.dot(inv_L_in,delta_z_in)

inv_L = np.linalg.inv(L)

delta_x = np.dot(inv_L,delta_z_in)

delta_z = np.dot(L,delta_x)

print('bins_to_remove:',bins_to_remove)
print('delta_z_in=',delta_z_in)
print('delta_x=',delta_x)
print('delta_x_in=',delta_x_in)
# print('delta_z=',delta_z)

nBins=5
for z in range(1,nBins+1):
	if z in bins_to_remove:
		print('removed bin:',z)
	else:
		print('delta_x_in-delta_x=','%0.3f' %  (delta_x_in-delta_x)[z-1])

print('\n\n')
print('bins_to_remove:',bins_to_remove)

print('')
for z in range(1,nBins+1):
	print('uncorr_bias_'+str(z),'=','-5.0','%0.3f' %delta_x[z-1],'5.0')



print('')
for z in range(1,nBins+1):
	print('uncorr_bias_'+str(z),'=','gaussian','%0.3f' %delta_x[z-1],'1.0')



