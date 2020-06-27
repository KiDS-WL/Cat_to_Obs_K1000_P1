import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.stats import chi2
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
from astropy.io import fits

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

#This makes episilons appear as epsilons rather than varepsilons
plt.rcParams["mathtext.fontset"] = "cm"

#set up the figure grid
tick_spacing = 2
fig,axes= plt.subplots(5,3,figsize=(8, 12),gridspec_kw={'hspace': 0, 'wspace': 0})

# Read in user input to set the patch, blind, zmin,zmax, nbootstrap
if len(sys.argv) <2: 
    print("Usage: %s LFVER LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid BLIND " % sys.argv[0])
    sys.exit(1)
else:
    LFVER=sys.argv[1] # catalogue version identifier
    BLIND=sys.argv[2] # catalogue version identifier

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
alphadata=ascii.read('%s/Calc_1pt_Stats/GeneralPlots/KAll.autocal.BlindA.alpha_VS_ZB.ZBcut0.1-1.2_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.txt'%MD)
alpha_mean=np.array(alphadata['alpha_1']+alphadata['alpha_2'])*0.5
alpha_err=np.sqrt(np.array(alphadata['err_alpha_1'])**2 + np.array(alphadata['err_alpha_2'])**2)*0.5 
alpha_low = alpha_mean-alpha_err
alpha_high = alpha_mean+alpha_err

Covfits='%s//data/kids/fits/xipm_KIDS1000_BlindA_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits'%MD

# read in the covariance
f = fits.open(Covfits)
COVMAT=f[1].data
nbins=9

# read in Csys data per tomo bin                                                                                                                                               
for iz in range(0,ntomobin):
    for jz in range(iz,ntomobin):

        gridpos=gridpos + 1
        
        #which grid cell do we want to plot this in?
        grid_x_E=np.int(gridpos/3)
        grid_y_E=gridpos % 3    #remainder
        ax=axes[grid_x_E,grid_y_E]

        tomochar='%s_%s'%(iz+1,jz+1)
        Csysfile='%s_%s_%s.dat'%(filetop,tomochar,LFVER)
        Csysdata = np.loadtxt(Csysfile)

        theta=Csysdata[:,1]
        Csys_p=Csysdata[:,3]
        Csys_m=Csysdata[:,4]
        err_Csys_p=Csysdata[:,5]
        err_Csys_m=Csysdata[:,6]
        epsfepsf=Csysdata[:,9]
        gepsf=Csysdata[:,12]

        npairs_weighted=Csysdata[:,20]

        #Lets see if we can't do something better for the CSys errors, as the Treecorr estimate doesn't
        #take the weights into account

        sigma_epsf=0.021667148635665007
        sigma_e=0.384   #both components

        err_gepsf=sigma_e*sigma_epsf/np.sqrt(npairs_weighted)
        err_epsfepsf=sigma_epsf*sigma_epsf/np.sqrt(npairs_weighted)

        an_err_Csys_p = Csys_p * np.sqrt(2*(err_gepsf/gepsf)**2 + (err_epsfepsf/epsfepsf)**2)
    
        # construct the model assuming constant alpha across the survey

        lilj=alpha_low[iz]*alpha_low[jz]
        lihj=alpha_low[iz]*alpha_high[jz]
        hilj=alpha_high[iz]*alpha_low[jz]
        hihj=alpha_high[iz]*alpha_high[jz]

        #Which is the minimum/maximum?
        mintestll_lh=np.minimum(lilj,lihj)
        mintesthl_hh=np.minimum(hilj,hihj)
        maxtesthh_lh=np.maximum(hihj,lihj)
        maxtesthl_hh=np.maximum(hilj,hihj)
        alal_min=np.minimum(mintestll_lh, mintesthl_hh)
        alal_max=np.maximum(maxtesthh_lh, maxtesthl_hh)
        
        xip_ss_low=alal_min*(Csysdata[:,9])
        xip_ss_high=alal_max*(Csysdata[:,9])

        psf_despf_low=2*alal_min*Csysdata[:,16]
        psf_despf_high=2*alal_min*Csysdata[:,16]

        depsf_depsf=Csysdata[:,17]

        tot_err_low =  xip_ss_low - alpha_high[iz]*Csysdata[:,16] - alpha_high[jz]*Csysdata[:,16] + depsf_depsf
        tot_err_high = xip_ss_high - alpha_low[iz]*Csysdata[:,16] - alpha_low[jz]*Csysdata[:,16] + depsf_depsf

        #read in the expectation value for the cosmic shear signal
        #xipdata=np.loadtxt('%s/ForBG/outputs/test_output_S8_fid_test/shear_xi_plus/bin_%d_%d.txt'%(MD,iz+1,ij+1))
        #thetatheory=np.loadtxt('%s/ForBG/outputs/test_output_S8_fid_test/shear_xi_plus/theta.txt'%(MD))*360.0*60.0/(2.0*np.pi)
        xiptheory=np.loadtxt('%s/ForBG/new_outputs/test_output_S8_fid_test/chain/output_test_A/shear_xi_plus_binned/bin_%d_%d.txt'%(MD,jz+1,iz+1))
        thetatheory=np.loadtxt('%s/ForBG/new_outputs/test_output_S8_fid_test/chain/output_test_A/shear_xi_plus_binned/theta_bin_%d_%d.txt'%(MD,jz+1,iz+1))
    
        # and plot the results with annotations of the bin combination and p-value
        ax.errorbar(theta, Csys_p/xiptheory, yerr=an_err_Csys_p/xiptheory, color='magenta',label='$\\xi_+^{\\rm sys}$')

        #ax.plot(theta, depsf_depsf/xiptheory, color='blue',label='$<\\delta \\epsilon^{PSF} \delta \\epsilon^{PSF}>$')

        #ax.fill_between(theta, xip_ss_low/xiptheory, xip_ss_high/xiptheory, color='lightgrey',label='$\\alpha_i \\alpha_j \, <\\epsilon^{PSF} \\epsilon^{PSF}>$',linestyle=':')
        
        #ax.fill_between(theta, tot_err_low/xiptheory, tot_err_high/xiptheory, color='lightgreen',label='$\\alpha_i \\alpha_j \, <\\epsilon^{PSF} \\epsilon^{PSF}> - (\\alpha_i + \\alpha_j)<\\epsilon^{PSF} \\delta \\epsilon^{PSF}> + <\\delta \\epsilon^{PSF} \delta \\epsilon^{PSF}>$',linestyle=':',alpha=0.5)
        
        #ax.plot(thetatheory,thetatheory*xipdata*0.03,color='blue',label='$3\%$ of $\\xi_+(\\theta)$')


        #we also want to include a band which shows the error
        istart=gridpos*nbins  #factor of 2 because of xi+ and xi-
        diag_for_izjz=np.diag(COVMAT[istart:istart+nbins,istart:istart+nbins])
        fact=0.1
        ax.fill_between(theta, fact*np.sqrt(diag_for_izjz)/xiptheory *-1.0, fact*np.sqrt(diag_for_izjz)/xiptheory, color='lightblue',label='$0.1 \\sigma_{\\xi_+}}$',linestyle=':',alpha=0.3)


        ax.annotate(tomochar, xy=(0.15,0.1),xycoords='axes fraction',
                    size=12, ha='right', va='top')
        #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        #ax.set_yscale('log')
        ax.set_xscale('log')        
        #ax.set_ylim(-9e-6,9e-6)
        if iz<1 or jz<2:
            ax.set_ylim(-0.3,0.3)
        elif iz==1:
            ax.set_ylim(-0.12,0.12)
        else:
            ax.set_ylim(-0.09,0.09)
        ax.set_xlim(0.5,300.0)
        ax.axhline(y=0, color='grey')
        ax.axhline(y=0.02, color='black', ls=':')
        ax.axhline(y=-0.02, color='black', ls=':')

        #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))  
        #ax.xaxis.set_major_locator(LogLocator(base=2))
        ax.set_xticks([1, 10, 100]) 
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        # only label the subplots at the edges
        ax.label_outer()

#add labels
axes[0,2].legend(fontsize=14,ncol=1,loc='upper right',frameon=True)
axes[2,0].set_ylabel('$\\xi_+^{\\rm sys}/\\xi_+^{\\Lambda {\\rm CDM}}$',fontsize=18)

axes[4,0].set_xlabel('$\\theta$ (arcmin)',fontsize=15)
axes[4,1].set_xlabel('$\\theta$ (arcmin)',fontsize=15)
axes[4,2].set_xlabel('$\\theta$ (arcmin)',fontsize=15)





outfile='figures/CSys_%s_auto_and_cross_BLIND_%s.png'%(LFVER,BLIND)
plt.tight_layout()
plt.savefig(outfile)
plt.show()


