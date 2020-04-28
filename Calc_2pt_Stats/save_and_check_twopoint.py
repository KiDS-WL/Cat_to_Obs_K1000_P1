

    ######################################
    ##  save_and_check_twopoint.py      ##
    ##  Chieh-An Lin                    ##
    ##  Version 2020.04.08              ##
    ######################################


import sys

import numpy as np
import scipy.interpolate as itp
import astropy.io.fits as fits

sys.path.append("kcap/modules/scale_cuts/")
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

###############################################################################
## Examples

def saveFitsTwoPoint_NOfZ_buceros():
    prefix = '/disk05/calin/91_Data/mockFootprint/buceros/MFP_selection/'
    nOfZNameList = [
        '%snOfZ_hist_BOSS_2dFLenS_combined_tomo0.dat' % prefix, 
        '%snOfZ_hist_BOSS_2dFLenS_combined_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo0.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo2.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo3.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo4.dat' % prefix
    ]
    nGalList     = [0.010839, 0.011899,   0.856520, 1.562567, 2.235793, 1.522305, 1.384634]
    sigmaEpsList = [0.284965, 0.278972, 0.288525, 0.279802, 0.289495]
    saveName = 'twoPoint_None_mean_None_cov_None_nOfZ_buceros.fits'
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_COSEBI_KV450():
    scDict = {'use_stats': 'En'.lower()}
    meanName = '/disk05/calin/91_Data/KiDS/KV-450_COSEBIs/En_nBins_5_0.50-300.00.ascii'
    cov = wtp2.makeCOSEBICov1()
    nOfZNameList, nGalList, sigmaEpsList = wtp2.makeKV450NOfZ()
    saveName = 'twoPoint_En_mean_KV450_cov_KV450_nOfZ_KV450.fits'
    
    saveFitsTwoPoint(
        nbTomoN=0, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=20,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='file', meanName=meanName,
        covTag='variable', covName=cov,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_copySameMean():
    scDict = {'use_stats': 'PeeE'.lower()}
    loadPrefix = '/disk05/calin/91_Data/mockFootprint/'
    savePrefix = ''
    
    ## Here, we grab the mean from an existing file
    mean = wtp.TwoPointWrapper.from_fits('/disk05/calin/91_Data/KiDS/kcap/Input_mean_cov_nOfZ/twoPoint_PeeE_mean_theoryNoise1_cov_theoryVultur_nOfZ_milvus.fits', covmat_name=None).makeMeanVector()
    
    covNameList  = [
        '%szosterops/MFP_for_others/cov_th_PeeE_obs.dat' % loadPrefix,
        '%svultur/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix,
        '%szosterops/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix
    ]
    nOfZNameList, nGalList, sigmaEpsList = wtp2.makeMilvusNOfZ(sourceOnly=True)
    saveNameList = [
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_theoryZosterops_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_simVultur_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_simZosterops_nOfZ_milvus.fits' % savePrefix
    ]
    
    for covName, saveName in zip(covNameList, saveNameList):
        saveFitsTwoPoint(
            nbTomoN=0, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS=None,
            scDict=scDict,
            meanTag='variable', meanName=mean,
            covTag='file', covName=covName,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

def saveFitsTwoPoint_list_withCut():
    scDict = {
        'use_stats': 'PneE PeeE'.lower(),
        'keep_ang_PneE'.lower(): '90 500',
        'keep_ang_PeeE'.lower(): '90 500'
    }
    covName   = '/disk05/calin/91_Data/mockFootprint/buceros/cosmokids/thps_cov_kids1000_buceros_bp_apo_simple_paircount_obs_bandpower_E_noap_list.dat'
    nOfZNameList, nGalList, sigmaEpsList = wtp2.makeBucerosNOfZ()
    saveName = 'twoPoint_BP_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_buceros.fits'
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS='/disk05/calin/91_Data/KiDS/kcap/Flinc/test_buceros/',
        scDict=scDict,
        meanTag='CosmoSIS', meanName=None,
        covTag='list', covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_Flinc():
    scDict = {
        'use_stats': 'gT gX xiP xiM PneE PneB PeeE PeeB En Bn'.lower()
    }
    nOfZNameList, nGalList, sigmaEpsList = wtp2.makeBucerosNOfZ()
    saveName = 'twoPoint_all_mean_theoryNoiseFree_cov_simBuceros_nOfZ_buceros.fits'
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc='/disk05/calin/91_Data/mockFootprint/',
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='Flinc', meanName='buceros',
        covTag='Flinc', covName='buceros',
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

###############################################################################

if __name__ == '__main__':
    saveFitsTwoPoint_NOfZ_buceros()
    print()
    saveFitsTwoPoint_COSEBI_KV450()
    print()
    saveFitsTwoPoint_copySameMean()
    print()
    saveFitsTwoPoint_list_withCut()
    
    ## This one won't work if you don't have the Flinc directory tree.
    #saveFitsTwoPoint_Flinc()

###############################################################################

