

    ######################################
    ##  save_and_check_twopoint.py      ##
    ##  Chieh-An Lin                    ##
    ##  Version 2020.04.08              ##
    ######################################


import sys

import numpy as np
import scipy.interpolate as itp
import astropy.io.fits as fits

sys.path.append("../../kcap/modules/scale_cuts/")
#import twopoint
import wrapper_twopoint as wtp
import wrapper_twopoint2 as wtp2


"""
Frequently Asked Questions

Q: Can I put only n(z) in?
A: Yes you can. Check `saveFitsTwoPoint_NOfZ_buceros`.

Q: Can I put only mean in?
A: Yes you can.

Q: Can I put only covariance in?
A: I don't know what for, but yes you can.

Q: I don't have n_gal or sigma_eps. Does it matter?
A: For chains, no. Just put dummy values.

Q: Can n(z) in a twopoint file be used by KCAP?
A: YES!!! In your KCAP ini, add something like:
     [fits_nz_lens]
     nz_file = %(TWOPT_FILENAME)s
     data_sets = LENS
     
     [fits_nz_source]
     nz_file = %(TWOPT_FILENAME)s
     data_sets = SOURCE

Q: How to add noise to a mean vector?
A: Not here. Go to your KCAP ini. In [load_data_and_cut_scales], add:
     simulate = T ;; Save a file?
     simulate_with_noise = T
     number_of_simulations = 1
     mock_filename = %{TWOPT_FILE_TO_BE_SAVED}s
   And run test sampler.
   
Q: How to grab the mean of a file and put it into my new file?
A: Check `saveFitsTwoPoint_copySameMean`.

Q: How to apply scale cuts?
A: Check `saveFitsTwoPoint_list_withCut`.

Q: How to use a list-formatted covariance?
A: Check `saveFitsTwoPoint_list_withCut`.

Q: How to put multiple statistics together?
A: Check `saveFitsTwoPoint_Flinc`.

Q: How to check my saved file?
A: Check `printTwoPointHDU`, `printTwoPoint_fromFile`, and `printTwoPoint_fromObj`.

Q: How to compare two files?
A: Check `unitaryTest`.
"""

###############################################################################
## Main function

def saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc='/disk05/calin/91_Data/mockFootprint/',
        prefix_CosmoSIS='/disk05/calin/91_Data/KiDS/kcap/Flinc/test_buceros/',
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nOfZNameList=None, nGalList=None, sigmaEpsList=None,
        saveName=None
    ):
    """
    This is a general function to save twopoint file.
    
    Parameters
    ----------
    nbTomoN : int, optional
        Number of lens bins
    nbTomoG : int, optional
        Number of source bins
    N_theta : int, optional
        Number of theta bins
    theta_min : float, optional
        Lower limit of theta bins
    theta_max : float, optional
        Upper limit of theta bins
    N_ell : int, optional
        Number of ell bins
    ell_min : float, optional
        Lower limit of ell bins
    ell_max : float, optional
        Upper limit of ell bins
    nbModes : int, optional
        Number of COSEBIs modes
    prefix_Flinc : string, optional
        Prefix of Flinc input directory; only concerned if meanTag = 'Flinc'
    prefix_CosmoSIS : string, optional
        Prefix of CosmoSIS theory input directory
        Only concerned if meanTag = 'CosmoSIS'
    scDict : dict, optional
        Dictionary containing scale-cut arguments
        Same format as in kcap ini files
        All dictionary keys & values have to be in lower case
    meanTag : {None, 'Flinc', 'CosmoSIS', 'variable', 'file'}, optional
        Method of mean input. One of
        
        None
            No mean vector
        
        ``Flinc``
            Calculate Flinc means specified by `prefix_Flinc`
            `meanName` is then interpreted as the bird tag to be used
        
        ``CosmoSIS``
            Read theory outputs specified by `prefix_CosmoSIS`
            `meanName` is then ignored
        
        ``variable``
            Read `meanName` directly as a python object,
            supposed to be a list or an array
            Should already be ordered
          
        ``file``
            Read `meanName` as the path to a single column file
            If the file does not have '.npy' or '.fits' as extension,
            it is interpreted as an ASCII file.
    meanName : object, optional
        See `meanTag`
    covTag : {None, 'Flinc', 'list', 'variable', 'file'}, optional
        Method of covariance input. One of
        
        None
            No covariance
        
        ``Flinc``
            Calculate Flinc covariance specified by `prefix_Flinc`
            `covName` is then interpreted as the bird tag to be used
        
        ``list``
            Read theory covariance specified by `covName` as under Benjamin's 
            list format
            If the file has several terms (G, NG, SSC, etc.), it sums up all 
            terms automatically.
            All nan values are automatically replaced with 0
            Should not already apply scale cuts in the file
        
        ``variable``
            Read `covName` directly as a python object,
            supposed to be a squared 2D-array
            Should already be ordered
          
        ``file``
            Read `covName` as the path to a file containing a squared matrix
            If the file does not have '.npy' or '.fits' as extension,
            it is interpreted as an ASCII file.
    nOfZNameList : None or string list
        List of n(z) file names
        One file for each tomographic bin
        Has to be ASCII
        Should share the same z bins
        If None, no n(z) will be saved
    nGalList : float list
        List of n_gal
    sigmaEpsList : float list
        List of sigma_eps
        The length can be nbTomoG or nbTomoN+nbTomoG
    savename : string
        Path of the twopoint file to be saved
        
    Returns
    -------
    Nothing, but output a file
    """
    
    wtp2.saveFitsTwoPoint(
        nbTomoN=nbTomoN, nbTomoG=nbTomoG,
        N_theta=N_theta, theta_min=theta_min, theta_max=theta_max,
        N_ell=N_ell, ell_min=ell_min, ell_max=ell_max,
        nbModes=nbModes,
        prefix_Flinc=prefix_Flinc,
        prefix_CosmoSIS=prefix_CosmoSIS,
        scDict=scDict,
        meanTag=meanTag, meanName=meanName,
        covTag=covTag, covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

# Reads in from the list of input_files and puts them all into a long vector. 
# Make sure that the ordering is correct, col starts from 1
def make_2pt_vector(input_files,col=1):
    for rp in range(len(input_files)):
        file= open(input_files[rp])
        data=np.loadtxt(file,comments='#')
        if rp==0:
            data_all=data[col-1].copy()
        else:
            data_all=np.hstack(data_all,data[col-1])
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
### Making fits files for Phase-1 real data

# Folder and file names for nofZ, for the sources it will depend on the blind
blind = 'A'
cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'

FolderNameNofZ = ''
FolderNameData = '/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/'
nBins_lens = 2
lens1 = 
lens2 = 

nBins_source = 5
source1 = 
source2 = 
source3 = 
source4 = 
source5 = 

# number density of galaxies per arcmin^2
nGal_lens = [ , ] 
# read from file
filename =
nGal_source = np.loadtxt(filename)

nGal_all = nGal_lens + nGal_source

# read from file
filename = 
sigma_e  = np.loadtxt(filename)

# read in the BP results
name = FolderNameData +'/Pgk/xi2bandpow_output_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_8_Ell_100.0_1500.0_zbins'
input_files=[]
cols = 2
for bin1 in range(nBins_lens):
    for bin2 in range(bin1,nBins_srouce):
        fileNameInput=FolderName+name+str(bin1+1)+str(bin2+1)+'.asc'
        input_files.append(fileNameInput)

name = FolderNameData +'/Pkk/xi2bandpow_output_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_8_Ell_100.0_1500.0_zbins'
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_srouce):
        fileNameInput=FolderName+name+str(bin1+1)+'_'+str(bin2+1)+'.dat'
        input_files.append(fileNameInput)

BP_vector_no_m_bias = make_2pt_vector(input_files,col=col)


# read in the COSEBIs results
name = FolderNameData+'/COSEBIS/En_COSEBIS_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_theta_0.5_300_zbins_'
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_srouce):
        fileNameInput=FolderName+name+str(bin1+1)+'_'+str(bin2+1)+'.asc'
        input_files.append(fileNameInput)

COSEBIs_vector_no_m_bias = make_2pt_vector(input_files)

# read in the xipm results, bin and sort
theta_min=0.5
theta_max=300.0
nTheta=9
name = FolderNameData+'/XI/XI_K1000_ALL_BLIND_'+blind+'_'+cat_version+'_nbins_4000_theta_0.5_300.0_zbins_'
for bin1 in range(nBins_source):
    for bin2 in range(bin1,nBins_srouce):
        fileNameInput=FolderName+name+str(bin1+1)+'_'+str(bin2+1)+'.asc'
        file= open(fileNameInput)
        xipm_in=np.loadtxt(file,comments='#')
        theta = xipm_in[:,0]
        xip   = xipm_in[:,3]
        xim   = xipm_in[:,4]
        weight= xipm_in[:,-1]
        xip_binned = rebin(theta,xip,weight,theta_min,theta_max,nTheta)
        xim_binned = rebin(theta,xim,weight,theta_min,theta_max,nTheta)
        if counter==1:
            xip_all = xip_binned.copy()
            xim_all = xim_binned.copy()
        else:
            xip_all = np.hstack(xip_all,xip_binned)
            xim_all = np.hstack(xim_all,xim_binned)

xipm_all = np.hstack(xip_all,xim_all)

# Make the data and Cov and redshift file for BP for KiDS1000 Phase-1
def saveFitsBP_list_KIDS1000():
    scDict = {
        'use_stats': 'PneE PeeE'.lower()
    }
    # where are the redshifts saved? They should have the same z.
    FolderName = FolderNameNofZ
    
    nOfZNameList = [ lens1 ,
                     lens2 ,
                     source1,
                     source2,
                     source3,
                     source4,
                     source5]

    nGalList     = nGal_all
    sigmaEpsList = sigma_e

    covName   = 'covariance/kids1000/BP_cov_list.dat'
    name_tag  = '_no_m_bias'
    saveName  = 'BP_KIDS1000_Blind'+blind+name_tag+'_'+cat_version+'.fits'
    
    saveFitsTwoPoint(
        nbTomoN=nBins_lens, nbTomoG=nBins_source,
        # N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        # nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='variable', meanName=BP_vector_no_m_bias,
        covTag='list', covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

###############################################################################
## Checks

def printTwoPointHDU(name, ind=1):
    """
    Print the content of a given HDU of a twopoint file
    """
    hdr  = fits.getheader(name, ind)
    data = fits.getdata(name, ind)
    
    print()
    print(hdr.tostring(sep='\n'))
    print(data)
    return

def printTwoPoint_fromFile(name):
    """
    Print the summary info of a twopoint file
    """
    HDUList = fits.open(name)
    
    ## Check default HDU
    print()
    print('Check default HDU:')
    
    if 'SIMPLE' in HDUList[0].header:
        HDUList = HDUList[1:]
        print('  Passed.')
    else:
        print('  No default HDU.')
        print('  Means that this file was not generated in the standard way.')
        print('  Will continue.')
    
    hdrList  = [HDU.header for HDU in HDUList]
    dataList = [HDU.data for HDU in HDUList]
    
    print()
    wtp2._checkExtensions_fromFile(hdrList)
    print()
    wtp2._checkCovariance_fromFile(hdrList)
    print()
    wtp2._checkSpectra_fromFile(hdrList)
    print()
    wtp2._checkKernels_fromFile(hdrList, dataList)
    print()
    wtp2._checkNGal_fromFile(hdrList)
    return 

def printTwoPoint(TP, mean=True, cov=True, nOfZ=True):
    """
    Print the summary info of a twopoint object
    Useful when you want to see what is really stocked in the python object
    """
    if mean:
        print()
        print('Spectra:')
        for spectrum in TP.spectra:
            print()
            wtp2._printSpectrum(spectrum)
    
    if cov:
        print()
        print('Covariance:')
        if hasattr(TP, 'covmat_info') and TP.covmat_info is not None:
            print()
            wtp2._printCovMatInfo(TP.covmat_info)
            print('Direct cov.shape = %s' % str(TP.covmat.shape))
        else:
            print()
            print('Did not find `covmat_info` attribute')
    
    if nOfZ:
        print()
        print('Kernels:')
        for kernel in TP.kernels:
            print()
            wtp2._printKernel(kernel)
    
    ##print(TP.windows)
    ##print(TP._spectrum_index)
    return

def printTwoPoint_fromObj(name, mean=True, cov=True, nOfZ=True):
    """
    Print the summary info of a twopoint file by reading it first as an object
    """
    try:
        TP = wtp.TwoPointWrapper.from_fits(name, covmat_name='COVMAT')
    except:
        TP = wtp.TwoPointWrapper.from_fits(name, covmat_name=None)
    printTwoPoint(TP, mean=mean, cov=cov, nOfZ=nOfZ)
    return

def unitaryTest(name1, name2):
    """
    Check if two files are strictly identical
    """
    wtp2.unitaryTest(name1, name2)
    return


