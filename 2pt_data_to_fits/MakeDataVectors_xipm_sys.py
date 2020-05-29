import numpy as np

# Reads in from the list of input_files and puts them all into a long vector. 
# Make sure that the ordering is correct, col starts from 1 instead of 0
def make_2pt_vector(input_files, m_corr,col=1):
    for rp in range(len(input_files)):
        file= open(input_files[rp])
        data=np.loadtxt(file,comments='#')
        if data.ndim==1:
            if rp==0:
                data_all      = data.copy()
                data_all_corr = data/m_corr[rp]
            else:
                data_all      = np.hstack((data_all,data))
                data_all_corr = np.hstack((data_all_corr,data/m_corr[rp]))
        else:
            if rp==0:
                data_all      = data[:,col-1].copy()
                data_all_corr = data[:,col-1]/m_corr[rp]
            else:
                data_all      = np.hstack((data_all,data[:,col-1]))
                data_all_corr = np.hstack((data_all_corr,data[:,col-1]/m_corr[rp]))
    return data_all,data_all_corr

def rebin(x,signal,weight,x_min,x_max,nbins):
    # print('rebinning now')
    binned_output=np.zeros((nbins,3))
    for ibins in range(nbins):
        x_binned=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+0.5))
        upperEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins+1.0))
        lowerEdge=np.exp(np.log(x_min)+np.log(x_max/x_min)/(nbins)*(ibins))
        good=((x<upperEdge) & (x>lowerEdge))
        # print(x_binned)
        if(good.any()):
            weight_sum=weight[good].sum()
            x_binned_weighted=(x[good]*weight[good]).sum()/weight_sum
            binned_output[ibins,0]=x_binned
            binned_output[ibins,1]=x_binned_weighted
            binned_output[ibins,2]=(signal[good]*weight[good]).sum()/weight_sum
            # print(ibins,weight_sum,len(weight[good]))
        else:
            print("WARNING: not enough bins to rebin to "+str(nbins)+" log bins")
    return binned_output


##################################################################################
### Making data vectors for Phase-1 real data

blind = 'A'

cat_shear_version = 'LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'
cat_mask_version  = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask'
cat_version = cat_mask_version+'_'+cat_shear_version

# This is were the raw data is saved on cuillin
FolderNameData = '/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/'
outputFolder   = "../data/kids/psf_systematic_corrected/"
nBins_lens     = 2
nBins_source   = 5

# fiducial values
filename="../data/kids/multiplicative_bias/Summary_multiplicative_Fid_unblinded.npy"
m=np.load(filename)[:,1]


#####################################################################################################
# XIPM
theta_min=0.5
theta_max=300.0
str_tmin='0.5'
str_tmax='300'
nTheta=9
counter=1
name_sys = FolderNameData+'CSys/CSys_5Z_'
name = FolderNameData+'/XI/XI_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_4000_theta_0.5_300.0_zbins_'
# xip_all_list =[]
# xim_all_list =[]
# xip_all_corr_list =[]
# xim_all_corr_list =[]
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_source):
        m_corr= (1.+m[bin2])*(1.+m[bin1])
        fileNameInput=name+str(bin1+1)+'_'+str(bin2+1)+'.asc'
        file= open(fileNameInput)
        xipm_in=np.loadtxt(file,comments='#')
        theta = xipm_in[:,0]
        xip   = xipm_in[:,3]
        xim   = xipm_in[:,4]
        weight= xipm_in[:,-1]
        xip_binned = rebin(theta,xip,weight,theta_min,theta_max,nTheta)
        xim_binned = rebin(theta,xim,weight,theta_min,theta_max,nTheta)
        # now read the sys files, these are already binned into 9 bins
        fileName_sys=name_sys+str(bin1+1)+'_'+str(bin2+1)+'_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.dat'
        file= open(fileName_sys)
        xipm_sys=np.loadtxt(file,comments='#')
        xip_sys   = xipm_sys[:,3]
        xim_sys   = xipm_sys[:,4]
        if(((bin1==3) & (bin2==3)) or ((bin1==3) & (bin2==4)) or ((bin1==4) & (bin2==4))):
            print(bin1+1,bin2+1)
            print(xip_sys/xip_binned[:,-1])
            xip_binned[:,-1] -= xip_sys
        if counter==1:
            xip_all = xip_binned[:,-1].copy()
            xim_all = xim_binned[:,-1].copy()
            xip_all_corr = xip_binned[:,-1]/m_corr
            xim_all_corr = xim_binned[:,-1]/m_corr
        else:
            xip_all = np.hstack((xip_all,xip_binned[:,-1]))
            xim_all = np.hstack((xim_all,xim_binned[:,-1]))
            xip_all_corr = np.hstack((xip_all_corr,xip_binned[:,-1]/m_corr))
            xim_all_corr = np.hstack((xim_all_corr,xim_binned[:,-1]/m_corr))
        counter+=1

        # xip_all_list.append(xip_binned[:,-1])
        # xim_all_list.append(xim_binned[:,-1])
        # xip_all_corr_list.append((xip_binned[:,-1]/m_corr))
        # xim_all_corr = np.hstack((xim_binned[:,-1]/m_corr))

# xip_all = np.asarray(xip_all_list)
# xim_all = np.asarray(xim_all_list)
# xip_all_corr = np.asarray(xip_all_corr_list)
# xim_all_corr = np.asarray(xim_all_corr_list)

xipm_all      = np.hstack((xip_all,xim_all))
xipm_all_corr = np.hstack((xip_all_corr,xim_all_corr))

name_tag = 'no_m_bias'
savename = outputFolder+'xipm_K1000_ALL_BLIND_'+blind+'_'+name_tag+'_'+cat_version+'_nbins_'+str(nTheta)+'_theta_'+str_tmin+'_'+str_tmax+'.asc'
np.savetxt(savename,xipm_all)

name_tag = 'with_m_bias'
savename = outputFolder+'xipm_K1000_ALL_BLIND_'+blind+'_'+name_tag+'_'+cat_version+'_nbins_'+str(nTheta)+'_theta_'+str_tmin+'_'+str_tmax+'.asc'
np.savetxt(savename,xipm_all_corr)
