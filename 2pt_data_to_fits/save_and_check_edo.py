
    ######################################
    ##  save_and_check_Phase1.py        ##
    ##  Marika Asgari                   ##
    ##  Version 2020.04.21              ##
    ######################################

# This is based on Linc's save_and_check_twopoint. 
# It has been adapted to make .fits files for the Phase1 real data

import sys

import numpy as np
import scipy.interpolate as itp
import astropy.io.fits as fits
import os

# set the path to scale_cuts here
sys.path.append("/net/home/fohlen13/dvornik/random_git/kcap/modules/scale_cuts/")
sys.path.append("/net/home/fohlen13/dvornik/random_git/Cat_to_Obs_K1000_P1/Calc_2pt_Stats/")
#import twopoint
import wrapper_twopoint as wtp
import wrapper_twopoint2 as wtp2

###############################################################################
## Main function

def saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5, nbObs=1,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        nnAuto=True,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nobsTag=None, nobsName=None,
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
        nbTomoN=nbTomoN, nbTomoG=nbTomoG, nbObs=1,
        N_theta=N_theta, theta_min=theta_min, theta_max=theta_max,
        N_ell=N_ell, ell_min=ell_min, ell_max=ell_max,
        nbModes=nbModes,
        nnAuto=True,
        prefix_Flinc=prefix_Flinc,
        prefix_CosmoSIS=prefix_CosmoSIS,
        scDict=scDict,
        meanTag=meanTag, meanName=meanName,
        covTag=covTag, covName=covName,
        nobsTag=nobsTag, nobsName=nobsName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return


# copied from cosmosis
def load_histogram_form(ext):
    # Load the various z columns.
    # The cosmosis code is expecting something it can spline
    # so  we need to give it more points than this - we will
    # give it the intermediate z values (which just look like a step
    # function)
    zlow = ext.data['Z_LOW']
    zhigh = ext.data['Z_HIGH']
    # First bin.
    i = 1
    bin_name = 'BIN{0}'.format(i)
    nz = []
    z = ext.data['Z_MID']
    # Load the n(z) columns, bin1, bin2, ...
    while bin_name in ext.data.names:
        col = ext.data[bin_name]
        nz.append(col)
        i += 1
        bin_name = 'BIN{0}'.format(i)
    # First bin.
    i = 1
    ngal_name = "NGAL_"+str(i)
    n_bar= []
    while ngal_name in ext.header.keys():
        n_b = ext.header[ngal_name]
        n_bar.append(n_b)
        i += 1
        ngal_name = "NGAL_"+str(i)
    nbin = len(nz)
    print("        Found {0} bins".format(nbin))
    nz = np.array(nz)
    # z, nz = ensure_starts_at_zero(z, nz)
    for col in nz:
        norm = np.trapz(col, z)
        col /= norm
    return z, nz, n_bar


def mkdir_mine(dirName):
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")


##################################################################################
### Making fits files for Phase-1 real data


# Folder and file names for nofZ, for the sources it will depend on the blind
blind = 'C'
cat_version = 'edo'
name_tag    = 'no_m_bias' # with_m_bias # no_m_bias # bmodes


bp_filename      = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/Pkids_data_full_ibc2.dat'


nBins_lens = 2
lens1 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz1.dat'
lens2 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz2.dat'

nBins_source = 4
source1 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_z1_kids_binned.dat'
source2 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_z2_kids_binned.dat'
source3 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_z3_kids_binned.dat'
source4 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_z4_kids_binned.dat'

smf = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/bin_13_mlf_full.txt'
smf2 = ['/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/bin_13_mlf_small_only.txt',
        '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/bin_13_mlf_large_only.txt']
 

z_sources_bins = np.linspace(0.0, 0.6, 50)#np.loadtxt(source1)[:,0]
nz_lens1 = np.loadtxt(lens1)
nz_lens1_new = np.interp(z_sources_bins, nz_lens1[:,0], nz_lens1[:,1], left=0.0, right=0.0)
np.savetxt('/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz1_new.dat', np.array([z_sources_bins, nz_lens1_new]).T)


nz_lens2 = np.loadtxt(lens2)
nz_lens2_new = np.interp(z_sources_bins, nz_lens2[:,0], nz_lens2[:,1], left=0.0, right=0.0)
np.savetxt('/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz2_new.dat', np.array([z_sources_bins, nz_lens2_new]).T)


lens1 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz1_new.dat'
lens2 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz2_new.dat'

obs1 = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/nz_GAMAz1_new.dat'

# number density of galaxies per arcmin^2

# Lenses:
n_2dflens    = np.asarray([0.002890,0.003674])
n_boss       = np.asarray([0.014478,0.016597])
area_2dflens = 342.879925
area_boss    = 322.255634
area_total   = 852.845901656631
nGal_lens_average = n_2dflens*area_2dflens/area_total+n_boss*area_boss/area_total
nGal_lens = [ nGal_lens_average[0],nGal_lens_average[1] ]

# Sources:
# read from file
filename = '/net/home/fohlen13/dvornik/random_git/Cat_to_Obs_K1000_P1/data/kids/number_density/ngal_blindC.ascii'
nGal_source = np.loadtxt(filename)[:4]

nGal_all = nGal_lens + (nGal_source).tolist()

# read from file
filename = '/net/home/fohlen13/dvornik/random_git/Cat_to_Obs_K1000_P1/data/kids/ellipticity_dispersion/sigma_e_blindC.ascii'
sigma_e  = np.loadtxt(filename)[:4]

# Make the data and Cov and redshift file for BP for KiDS1000 Phase-1
def saveFitsBP_list_KIDS1000():
    scDict = {
        'use_stats': 'PeeE PneE Pnn 1pt'.lower()
    }
    
    nOfZNameList = [ lens1 ,
                     lens2 ,
                     source1,
                     source2,
                     source3,
                     source4,
                     obs1,]

    nGalList     = nGal_all
    sigmaEpsList = sigma_e.tolist()
    
    covName   = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/Pkidscov_data_PlanckIter2.dat'
    
    saveName  = '/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/data_edo_test_tmp.fits'
    
    saveFitsTwoPoint(
        nbTomoN=nBins_lens, nbTomoG=nBins_source, nbObs=1,
        N_ell=5, ell_min=155, ell_max=1930,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        nnAuto=True,
        meanTag='file', meanName=bp_filename,
        covTag='name', covName=covName,
        nobsTag='file', nobsName=smf,
        nOfZNameList=nOfZNameList,# nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return


###############################################################################
## Checks and plots

def plot_redshift(filename,title,savename):
    import matplotlib.pyplot as plt
    F=fits.open(filename)
    ext=F["nz_source"]
    z_source, nz_source, n_bar_source=load_histogram_form(ext)

    try:
        ext=F["nz_lens"]
        z_lens, nz_lens, n_bar_lens=load_histogram_form(ext)
        plot_lenses=True
    except:
        print("no lenses given")
        plot_lenses=False

    F.close()

    if(plot_lenses):
        plt.clf()
        ax=plt.subplot(2,1,1)
        plt.ylabel("P(z)")
        plt.title(title)
        plt.setp(ax.get_xticklabels(),  visible=False)
        plt.subplots_adjust(wspace=0,hspace=0)
        for bin1 in range(len(nz_lens)):
            plt.xlim(0,2.0)
            #print(z_lens)
            plt.plot(z_lens,nz_lens[bin1],label='lens '+str(bin1+1))
            plt.legend(loc='best')

        ax=plt.subplot(2,1,2)
        plt.setp(ax.get_xticklabels(),  visible=True)
        for bin1 in range(len(nz_source)):
            plt.xlim(0,2.0)
            plt.plot(z_source,nz_source[bin1],label='source '+str(bin1+1))
            plt.legend(loc='best')

        plt.xlabel("z")
        plt.ylabel("P(z)")
        plt.savefig(savename,bbox_inches='tight')
    else:
        plt.clf()
        ax=plt.subplot(1,1,1)
        plt.setp(ax.get_xticklabels(),  visible=True)
        for bin1 in range(len(nz_source)):
            plt.xlim(0,2.0)
            plt.plot(z_source,nz_source[bin1],label='source '+str(bin1+1))
            plt.legend(loc='best')

        plt.xlabel("z")
        plt.ylabel("P(z)")
        plt.savefig(savename,bbox_inches='tight')


def plot_covariance(filename,title,savename):
    import matplotlib.pyplot as plt
    F=fits.open(filename)
    ext=F["COVMAT"]
    covariance= ext.data
    fig, ax = plt.subplots()
    im = ax.imshow(covariance)
    cbar = ax.figure.colorbar(im, ax=ax)
    plt.title(title)
    plt.savefig(savename)

def plot_correlation_mat(filename,title,savename):
    import matplotlib.pyplot as plt
    F=fits.open(filename)
    ext=F["COVMAT"]
    cov= ext.data
    corr=np.zeros((len(cov),len(cov)))
    for i in range(len(cov)):
        for j in range(len(cov)):
            corr[i,j]=cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])
    fig, ax = plt.subplots()
    im = ax.imshow(corr)
    cbar = ax.figure.colorbar(im, ax=ax)
    plt.title(title)
    plt.savefig(savename)

def plot_data(filename,title,extname,savename):
    import matplotlib.pyplot as plt
    F=fits.open(filename)
    ext=F[extname]
    data=ext.data['VALUE']
    x_index = ext.data['ANGBIN']
    x_val   = ext.data['ANG']
    plt.clf()
    plt.title(title)
    plt.plot(data,'x')
    plt.savefig(savename)
    
def plot_data_1pt(filename,title,extname,savename):
    import matplotlib.pyplot as plt
    F=fits.open(filename)
    ext=F[extname]
    #print(ext.columns)
    plt.clf()
    plt.title(title)
    for i in range(1, 99999):
        try:
            data = ext.data['BIN{}'.format(i)]
            x_val = ext.data['OBS{}'.format(i)]
            #print(x_val, data)
            plt.plot(x_val,data)
        except:
            break
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(savename)


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


############################################################
# plot things here

saveFitsBP_list_KIDS1000()

# exit()
FolderPlots='/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/plots'
mkdir_mine(FolderPlots)




# BP
filename='/net/home/fohlen13/dvornik/halo_model_mc/tests/pmm/test_run/vUitert18_KiDS450_data/data_edo_test_tmp.fits'
title= 'KiDS450-GAMA'
savename=FolderPlots+'/KiDS1000_nofz_test_tmp.pdf'
plot_redshift(filename,title,savename)

savename=FolderPlots+'/BP_covariance_test_tmp.pdf'
plot_covariance(filename,title,savename)

savename=FolderPlots+'/BP_correlation_matrix_test_tmp.pdf'
plot_correlation_mat(filename,title,savename)

#file=open(bp_filename)
#bp=np.loadtxt(file)

extname='PeeE'
savename=FolderPlots+'/BP_data'+extname+'_test_tmp.pdf'
plot_data(filename,title,extname,savename)

extname='PneE'
savename=FolderPlots+'/BP_data'+extname+'_test_tmp.pdf'
plot_data(filename,title,extname,savename)

extname='Pnn'
savename=FolderPlots+'/BP_data'+extname+'_test_tmp.pdf'
plot_data(filename,title,extname,savename)

extname='1pt'
savename=FolderPlots+'/BP_data'+extname+'_test_tmp.pdf'
plot_data_1pt(filename,title,extname,savename)



