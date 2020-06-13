import numpy as np

# Reads in from the list of input_files and puts them all into a long vector. 
# Make sure that the ordering is correct, col starts from 1 instead of 0
def make_2pt_vector(input_files,col=1):
    for rp in range(len(input_files)):
        file= open(input_files[rp])
        data=np.loadtxt(file,comments='#')
        if data.ndim==1:
            if rp==0:
                data_all      = data.copy()
            else:
                data_all      = np.hstack((data_all,data))
        else:
            if rp==0:
                data_all      = data[:,col-1].copy()
            else:
                data_all      = np.hstack((data_all,data[:,col-1]))
    return data_all

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

blind = 'C'
cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'

# This is were the raw data is saved on cuillin
FolderNameData = '/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/'
outputFolder   = "../data/kids/"
nBins_lens     = 2
nBins_source   = 5

#####################################################################################################
# BP
name = FolderNameData +'/Pgk/xi2bandpow_output_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_8_Ell_100.0_1500.0_zbins_'
input_files = []
col = 4
for bin1 in range(nBins_lens):
    for bin2 in range(nBins_source):
        fileNameInput=name+str(bin1+1)+'_'+str(bin2+1)+'.dat'
        input_files.append(fileNameInput)


name = FolderNameData +'/Pkk/xi2bandpow_output_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_8_Ell_100.0_1500.0_zbins_'
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_source):
        fileNameInput=name+str(bin1+1)+'_'+str(bin2+1)+'.dat'
        input_files.append(fileNameInput)

BP_vector = make_2pt_vector(input_files,col=col)

name_tag ='bmodes'

savename = outputFolder+'bp_K1000_ALL_BLIND_'+blind+'_'+name_tag+'_'+cat_version+'_nbins_8_Ell_100.0_1500.0.asc'
np.savetxt(savename,BP_vector)


#####################################################################################################
# COSEBIs
input_files = []
name = FolderNameData+'/COSEBIS/Bn_COSEBIS_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_theta_0.5_300_zbins_'
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_source):
        fileNameInput=name+str(bin1+1)+'_'+str(bin2+1)+'.asc'
        input_files.append(fileNameInput)

COSEBIs_vector  = make_2pt_vector(input_files)

name_tag = 'bmodes'
savename = outputFolder+'cosebis_K1000_ALL_BLIND_'+blind+'_'+name_tag+'_'+cat_version+'_nbins_theta_0.5_300.asc'
np.savetxt(savename,COSEBIs_vector_with_m_bias)

