# 07/04/2020, B. M. Giblin, Postdoc, Edinburgh
# Read in and compare the model fits which include and exclude
# the impact of mgnification

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.gridspec as gridspec

plt.rc('font', size=10)
# ----- Load Input ----- #                                                                                                               
from Get_Input import Get_Input
paramfile = sys.argv[1]   # e.g. params_KiDSdata.dat                                                                                     
GI = Get_Input(paramfile)
SOURCE_TYPE = GI.Source_Type()
LENS_TYPE = GI.Lens_Type()
RANDOM_TYPE = GI.Random_Type()

Cov_Method = "Spin"   # The method for calculating the gamma_t realisations for use in covariance estimation
INDIR='Output/SOURCE-%s_LENS-%s' %(SOURCE_TYPE, LENS_TYPE)
if "MICE2" in SOURCE_TYPE:
    # Additional identifiers when working with MICE
    # True/Estimated P(z) and Magnification on/off
    Mag_OnOff = GI.Mag_OnOff()
    Pz = GI.Pz_TrueEstimated()
    SN = GI.SN()
    INDIR += '_Pz%s_SN%s_mag%s' %(Pz,SN,Mag_OnOff)

OUTDIR = INDIR + '/' + Cov_Method.upper()

params_m = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True))
params_nm = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False))
thetas, model_m = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True),
                           usecols=(0,1), unpack=True)
thetas, model_nm = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False),
                           usecols=(0,1), unpack=True)

nz_source = 5
nz_lens = 5
nz_bins = nz_source*nz_lens
ntheta = 4

# Plot the model
def Plot_gt_MagOnOff(yaxis_label):
    fig = plt.figure(figsize = (12,10)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, nz_lens)
    k=0
    for i in range(nz_source):
        for j in range(nz_lens):
            ax1 = plt.subplot(gs[k])
            plt.plot( thetas[k*ntheta:(k+1)*ntheta], thetas[k*ntheta:(k+1)*ntheta]*model_m[k*ntheta:(k+1)*ntheta],
                      color='red', label='Mag On' )
            plt.plot( thetas[k*ntheta:(k+1)*ntheta], thetas[k*ntheta:(k+1)*ntheta]*model_nm[k*ntheta:(k+1)*ntheta],
                      color='blue', label='Mag Off' )

            ax1.set_ylim([-0.005, 0.02])

            if i==nz_source-1:
                # Get an x-axis for bottom-row.
                ax1.set_xlabel(r'$\theta$ [arcmin]')
                if j==nz_lens-1:
                    # Get a legend just for the very bottom-right subplot.
                    ax1.legend(loc='upper left', frameon=False)
            else:
                ax1.set_xticks([])

            if j==0:
                # Get a y-axis for first column
                ax1.set_ylabel(yaxis_label)
            else:
                ax1.set_yticks([])
            k+=1
            ax1.text(0.95,0.95, 't%s, sp%s'%(i+1,j+1), ha="right", va="top", transform=plt.gca().transAxes)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    return
Plot_gt_MagOnOff(r'$\theta\gamma_{\rm t}^{\rm model}$')


def Plot_Amplitudes():
    plt.figure()
    #plt.scatter( np.arange(20)+1, 1e3*params_m[:20], marker='x', color='red', label='Mag On' )
    #plt.scatter( np.arange(20)+1, 1e3*params_nm[:20], marker='o', color='blue', label='Mag Off' )

    plt.scatter( np.arange(20)+1, params_m[:20]/params_nm[:20], marker='o', color='black')
    plt.xlabel( r'bin number' )
    #plt.ylabel( r'Amplitude $[\times 10^{-3}]$' )
    plt.ylabel( r'Amplitude ratio: mag/no-mag' )
    plt.title( r'Amplitude ratio: %.2f pm %.2f' %( np.mean(params_m[:20]/params_nm[:20]),
                                                   np.std( params_m[:20]/params_nm[:20])))
    plt.legend(loc='best', frameon=True)
    plt.show()
    return
Plot_Amplitudes()
    


        
