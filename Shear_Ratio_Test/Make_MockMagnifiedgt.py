# 07/04/2020, B. M. Giblin, Postdoc, Edinburgh
# Make a mock gamma_t data vector that has some level
# of magnification added in.


import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.gridspec as gridspec
import os

plt.rc('font', size=10)
# ----- Load Input ----- #

from Get_Input import Get_Input
paramfile = sys.argv[1]   

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

OUTDIR = INDIR.split('SOURCE-')[0] + 'SOURCE-TOY_' + INDIR.split('SOURCE-')[-1] # +Cov_Method.upper()
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# Make mock values of alpha (chosen randomly - but think mag should increase with lens bin)
alphas  = np.array([ 4.,4.,4.,4.,4. ]) #1.1, 1.3, 1.6, 1.9, 2.4 ])
Magnif_Shape = np.load('/home/bengib/kcap_NewInst/kcap/examples/output_magnification_alpha3.0/SRTparam_Bij.npy')

    
# Read in BOTH the no-mag measurement from the LENS/SOURCE specified,
# AND the no-mag model fit to this measurement (so we can chose which
# is the result we add magnification too and save as our toy model).
nz_source = 5
nz_lens = 5
nz_bins = nz_source*nz_lens
ntheta = 4

# These arrays are to store the no-mag/mag versions of the measurement
gt_nm = np.zeros([ nz_bins, ntheta ])
gt_m = np.zeros([ nz_bins, ntheta ])

# These arrays are the no-mag/mag versions of the model fit
gt_nm_mod = np.loadtxt('%s/SPIN/%sx%s_FitModel_MagFalse.dat'%(INDIR,SOURCE_TYPE,LENS_TYPE),
                                       usecols=(1,), unpack=True)
gt_nm_mod = np.reshape(gt_nm_mod, (nz_bins, ntheta))
gt_m_mod = np.zeros([ nz_bins, ntheta ])

k=0
for i in range(nz_source):
    for j in range(nz_lens):
        theta, gt_nm[k,:] = np.loadtxt('%s/GT_6Z_source_%s_5Z_lens_%s.asc' %(INDIR,i+1,j+1),
                                    usecols=(1,3), unpack=True)
        gt_m[k,:] = gt_nm[k,:] + 2*(alphas[j]-1.)*Magnif_Shape[k,:]
        gt_m_mod[k,:] = gt_nm_mod[k,:] + 2*(alphas[j]-1.)*Magnif_Shape[k,:]
        # Save the gt_m mock array:
        np.savetxt('%s/GT_6Z_source_%s_5Z_lens_%s.asc' %(OUTDIR,i+1,j+1),
                   np.c_[theta,gt_m_mod[k,:], gt_m[k,:],gt_m[k,:]], header='meanr gamT gamX sigma')
        k+=1
        
        
# Plot the gamma_t
def Plot_gt_MagOnOff(yaxis_label):
    fig = plt.figure(figsize = (12,10)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, nz_lens)
    k=0
    for i in range(nz_source):
        for j in range(nz_lens):
            ax1 = plt.subplot(gs[k])
            plt.plot( theta, theta*gt_nm_mod[k,:], color='blue', label='Mag Off' )
            plt.plot( theta, theta*gt_m_mod[k,:], color='red', label='Mag On' )

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
Plot_gt_MagOnOff(r'$\theta\gamma_{\rm t}$')
