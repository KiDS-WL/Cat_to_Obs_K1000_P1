# 07/04/2020, B. M. Giblin, Postdoc, Edinburgh
# Read in and compare the model fits:

Compare_Mag = False         # Compare fit params/model including/excluding
                           # the impact of magnification.
Compare_Single_Bin = True  # Compare fit params obtained for each l+s bin individually
                           # across lens bins (should be consistent)

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

nz_source = 5
nz_lens = 5
nz_bins = nz_source*nz_lens
ntheta = 4


if Compare_Mag:
    # Compare the fit params and model for Mag On with Mag Off
    params_m = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True))
    params_nm = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False))
    thetas, model_m = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True),
                                 usecols=(0,1), unpack=True)
    thetas, model_nm = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False),
                                  usecols=(0,1), unpack=True)


    
if Compare_Single_Bin:
    # Compare the fit params obtained for each l+s bin, across lens bins
    params_sl = np.zeros([nz_source, nz_lens])
    params_sl_err = np.zeros([nz_source, nz_lens])
    for i in range(nz_source):
        for j in range(nz_lens):
            params_sl[i,j] = np.loadtxt(OUTDIR+'/%sx%s_FitParams_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j))
            params_sl_err[i,j] = 2*np.loadtxt(OUTDIR+'/%sx%s_FitParamsErr_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j))
            
    #thetas, model_m = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True),
    #                             usecols=(0,1), unpack=True)

# ------------------------------------- FUNCTIONS FOR COMPARE MAG ON/OFF ANALYSES -------------------------------------

def Plot_gt_MagOnOff(yaxis_label):
    # Plot the Mag On/Off best-fit models
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

def Plot_MagOnOff_Amplitudes():
    # Compare the amplitudes obtained in the Mag On/Off cases
    
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

if Compare_Mag:
    Plot_gt_MagOnOff(r'$\theta\gamma_{\rm t}^{\rm model}$')
    Plot_MagOnOff_Amplitudes()
    
# ------------------------------------- ^^^FUNCTIONS FOR COMPARE MAG ON/OFF ANALYSES^^^ -------------------------------------





# ------------------------------------- FUNCTIONS FOR COMPARE SINGLE BIN ANALYSIS -------------------------------------
def Plot_SingleBin_Params():
    ylimits = [ [ -2.5, 3. ], [ 0.0, 0.08], [ 0.0, 0.059], [ 0.0, 0.039], [ 0.015, 0.039] ]
    
    # Plot the Mag On/Off best-fit models
    fig = plt.figure(figsize = (8,6)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, 1)
    colors=[ 'red', 'blue', 'orange', 'green', 'magenta' ]
    for i in range(nz_source):
        ax1 = plt.subplot(gs[i])
        ax1.fill_between( np.arange(6)+0.5 , (params_sl[i,0]-params_sl_err[i,0]),
                          (params_sl[i,0]+params_sl_err[i,0]), color='dimgrey', alpha=0.5 )
        ax1.errorbar( range(1,nz_lens+1), params_sl[i,:], yerr=params_sl_err[i,:], fmt='o',
                      color=colors[i], edgecolor=None )
        if i==(nz_source-1):
            ax1.set_xlabel('Tomo bin number')
            #ax1.set_xticks([])
        else:
            ax1.set_xticks([])

        if i==2:
            ax1.set_ylabel('Amplitude')
        
        ax1.set_ylim([ ylimits[i][0] , ylimits[i][1] ])
        #ax1.set_ylim([ 0.5*(params_sl[i,0]-params_sl_err[i,0]).min(),
        #               2.0*(params_sl[i,0]+params_sl_err[i,0]).max() ])
        ax1.text(0.95,0.95, 'sp%s'%(i+1), ha="right", va="top", transform=plt.gca().transAxes)
        
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    return
Plot_SingleBin_Params()
# ------------------------------------- ^^^FUNCTIONS FOR COMPARE SINGLE BIN ANALYSIS^^^ ------------------------------------- 
