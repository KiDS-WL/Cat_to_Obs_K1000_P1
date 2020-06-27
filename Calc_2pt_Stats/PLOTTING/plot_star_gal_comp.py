import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.stats import chi2
import matplotlib.ticker as ticker
from astropy.io import fits

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)

#This makes episilons appear as epsilons rather than varepsilons
plt.rcParams["mathtext.fontset"] = "cm"

#set up the figure grid
tick_spacing = 2
fig,axes= plt.subplots(5,1,figsize=(6.5, 10),gridspec_kw={'hspace': 0, 'wspace': 0})

# Read in user input to set the patch, blind, zmin,zmax, nbootstrap
if len(sys.argv) <2: 
    print("Usage: %s LFVER LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid BLIND" % sys.argv[0])
    sys.exit(1)
else:
    LFVER=sys.argv[1] # catalogue version identifier
    BLIND=sys.argv[2] # blind

# number of tomographic bins
ntomobin=5

# tomake different sub plots we count the grid square that we want to plot in
#initialising the counter
gridpos=-1

#information about the file names
#filetop='/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/COSEBIS/Bn_COSEBIS_K1000_ALL_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c'
#filetop='/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/CSys/CSys_5Z'
filetop='CSyS/CSys_BLIND_%s_5Z'%BLIND
#_1_LF_svn_309c_2Dbins.dat

#MD='/home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/'
#MD='/Users/heymans/KiDS/Cat_to_Obs_K1000_P1/'

MD='/Users/macleod/CH_work/Cat_to_Obs_K1000_P1/'
# Read in alphas from 1pt Output data directory
#alphadata=ascii.read('%s/Calc_1pt_Stats/Output_data/KAll.autocal.BlindA.alpha_VS_ZB.ZBcut0.1-1.2_%s.dat'%(MD,LFVER))
#alphadata=ascii.read('%s/Calc_1pt_Stats/GeneralPlots/KAll.autocal.BlindA.alpha_VS_ZB.ZBcut0.1-1.2_%s_THELI_INT.txt'%(MD,LFVER))
alphadata=ascii.read('%s/Calc_1pt_Stats/GeneralPlots/KAll.autocal.Blind%s.alpha_VS_ZB.ZBcut0.1-1.2_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.txt'%(MD,BLIND))
alpha_mean=np.array(alphadata['alpha_1']+alphadata['alpha_2'])*0.5
alpha_err=np.sqrt(np.array(alphadata['err_alpha_1'])**2 + np.array(alphadata['err_alpha_2'])**2)*0.5 
alpha_low = alpha_mean-2*alpha_err
alpha_high = alpha_mean+2*alpha_err

Covfits='%s//data/kids/fits/xipm_KIDS1000_BlindA_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits'%MD

# read in the covariance
f = fits.open(Covfits)
COVMAT=f[1].data
nbins=9

# read in Csys data per tomo bin                                                                                                                                               
for iz in range(0,ntomobin):
    gridpos=gridpos + 1
        
    #which grid cell do we want to plot this in?
    grid_x_E=np.int(gridpos/3)
    grid_y_E=gridpos % 3    #remainder
    ax=axes[iz]

    tomochar='%s_%s'%(iz+1,iz+1)
    Csysfile='%s_%s_%s.dat'%(filetop,tomochar,LFVER)
    Csysdata = np.loadtxt(Csysfile)

    theta=Csysdata[:,1]
    Csys_p=Csysdata[:,3]
    Csys_m=Csysdata[:,4]
    err_Csys_p=Csysdata[:,5]
    err_Csys_m=Csysdata[:,6]

    epsf_epsf=Csysdata[:,9]
    g_epsf=Csysdata[:,12] 
    epsf_depsf=Csysdata[:,16] 
    depsf_depsf=Csysdata[:,17]   
    g_depsf=Csysdata[:,18]   
    err_g_depsf=Csysdata[:,19] 
    npairs_weighted=Csysdata[:,20]

    # can we do an analytical estimate of the errors
    # the error on delta_epsf is zero
    # so the error will be
    sigma_depsf=0.0007459803415640388
    sigma_epsf=0.021667148635665007
    sigma_e=0.384   #both components
    average_epsf_sq=5.0310059452828268e-05
    gdepsf_analytical_error = sigma_e*sigma_depsf/np.sqrt(npairs_weighted)
    gepsf_analytical_error = sigma_e*sigma_epsf/np.sqrt(npairs_weighted)

    g_epsf_mod_low=alpha_low[iz]*epsf_epsf -alpha_high[iz]*average_epsf_sq  - epsf_depsf 
    g_epsf_mod_high=alpha_high[iz]*epsf_epsf -alpha_low[iz]*average_epsf_sq - epsf_depsf

    g_depsf_mod_low=(alpha_low[iz]*epsf_depsf - depsf_depsf)*10 
    g_depsf_mod_high=(alpha_high[iz]*epsf_depsf - depsf_depsf)*10 

    # and plot the results with annotations of the bin combination and p-value
    #ax.errorbar(theta, Csys_p/xiptheory, yerr=err_Csys_p/xiptheory, color='magenta',label='$\\xi_+^{sys}$')

    MF=1e5
    
    ax.fill_between(theta, theta**0.5*g_epsf_mod_low*MF, theta**0.5*g_epsf_mod_high*MF, color='lightblue',label='${\\rm Eq. 15}$',linestyle=':')
    #label='$\\alpha \\langle\\epsilon^{PSF} \\epsilon^{PSF}\\rangle  - \\langle\\epsilon^{PSF} \\delta\\epsilon^{PSF}\\rangle - \\alpha [\\overline{\\epsilon^{PSF}}]^2$',linestyle=':')

    ax.fill_between(theta, theta**0.5*g_depsf_mod_low*MF, theta**0.5*g_depsf_mod_high*MF, color='grey',label='$ 10 \\times \, {\\rm Eq. 16}$',linestyle=':')
    #label='$ 10\,\\left(\\alpha \\langle\\epsilon^{PSF} \\delta\\epsilon^{PSF}\\rangle - \\langle\\delta\\epsilon^{PSF} \\delta \\epsilon^{PSF}\\rangle - \\alpha \\overline{\\epsilon^{PSF}} \, \\overline{\\delta\\epsilon^{PSF}}\,\\right)$'

    ax.errorbar(theta, theta**0.5*g_epsf*MF, yerr=theta**0.5*gepsf_analytical_error*MF, color='blue',label='$\\langle\\epsilon^{\\rm obs} \\epsilon^{\\rm PSF}\\rangle$')
    ax.errorbar(theta, theta**0.5*g_depsf*10*MF, yerr=theta**0.5*gdepsf_analytical_error*10*MF, color='magenta',label='$10 \\langle\\epsilon^{\\rm obs} \\delta \\epsilon^{\\rm PSF}\\rangle$')
    
    #    ax.fill_between(theta, tot_err_low/xiptheory, tot_err_high/xiptheory, color='lightgreen',label='toterr',linestyle=':',alpha=0.5)
        
        #ax.plot(thetatheory,thetatheory*xipdata*0.03,color='blue',label='$3\%$ of $\\xi_+(\\theta)$')

    ax.annotate(tomochar, xy=(0.09,0.15),xycoords='axes fraction',
                    size=12, ha='right', va='top')
        #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        #ax.set_yscale('log')
    ax.set_xscale('log')        
    ax.set_ylim(-4.9e-5*MF,4.9e-5*MF)

    #if iz<1 :
    #    ax.set_ylim(-0.3,0.3)
    #elif iz==1:
    #    ax.set_ylim(-0.12,0.12)
    #else:
    #    ax.set_ylim(-0.09,0.09)
    ax.set_xlim(0.5,300.0)
    ax.axhline(y=0, color='black', ls=':')
    ax.axhline(y=0.02, color='black', ls=':')
    ax.axhline(y=-0.02, color='black', ls=':')

    # only label the subplots at the edges
    ax.label_outer()

#add labels
axes[0].legend(fontsize=15,ncol=2,loc='upper center',frameon=False,bbox_to_anchor=(0.42, 1.6))
axes[2].set_ylabel('$\\xi_+ \\sqrt{\\theta}\,\,\,\,\,\,  [10^{-5} \, {\\rm arcmin}^{0.5}\,]$',fontsize=16)
axes[4].set_xlabel('$\\theta \,\,({\\rm arcmin})$',fontsize=16)






outfile='figures/star_gal_comp_%s_auto_and_cross_BLIND_%s.png'%(LFVER,BLIND)
plt.tight_layout()
plt.savefig(outfile)
plt.show()


