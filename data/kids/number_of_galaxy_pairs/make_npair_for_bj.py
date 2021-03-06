import numpy as np

def rebin_sum(x,signal,x_min,x_max,nbins):
	# print('rebinning now')
	binned_output=np.zeros((nbins,2))
	for ibins in range(nbins):
		x_binned=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+0.5))
		upperEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+1.0))
		lowerEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins))
		good=((x<upperEdge) & (x>lowerEdge))
		# print(x_binned)
		if(good.any()):
			binned_output[ibins,0]=x_binned
			binned_output[ibins,1]=(signal[good]).sum()
		else:
			print("WARNING: not enough bins to rebin to "+str(nbins)+" log bins")
	return binned_output

nTheta    = 9
theta_min = 0.5
theta_max = 300.0
nBins_lens   = 2
nBins_source = 5
n_r_pair_clustering = int(nBins_lens*(nBins_lens+1)/2.)
n_r_pair_ggl = int(nBins_lens*nBins_source)
n_r_cosmicshear= int(nBins_source*(nBins_source+1)/2.)
n_redshift_pairs= int(n_r_pair_clustering+n_r_pair_ggl+n_r_cosmicshear)
npair_mat = np.ones((nTheta,n_redshift_pairs))

cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'

blind = 'C'

# We are going to only put values for the cosmic shear columns, the rest are set to one

start_rp = n_r_pair_ggl+n_r_pair_clustering
index = 0 
for z1 in range(nBins_source):
	for z2 in range(z1, nBins_source):
		filename='../xipm/XI_K1000_ALL_BlIND_'+blind+'_'+cat_version+'_nbins_4000_theta_0.5_300.0'+'_nBins_5_Bin'+str(z1+1)+'_Bin'+str(z2+1)+'.ascii'
		xipm=np.loadtxt(filename,comments='#')
		theta_all = xipm[:,0]
		npair_all=xipm[:,-1]
		# now bin it
		npair_binned=rebin_sum(theta_all,npair_all,theta_min,theta_max,nTheta)
		npair_mat[:,start_rp+index] = npair_binned[:,1]
		index+=1

filename = 'npair_blind'+blind+'_'+cat_version+'_nTheta'+str(nTheta)+'.ascii'
np.savetxt(filename,npair_mat)


filename = 'theta_blind'+blind+'_'+cat_version+'_nTheta'+str(nTheta)+'.ascii'
np.savetxt(filename,npair_binned[:,0])




		