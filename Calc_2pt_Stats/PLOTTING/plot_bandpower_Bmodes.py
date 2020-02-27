import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.stats import chi2
import matplotlib.ticker as ticker


# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}

plt.rc('font', **font)

#set up the figure grid
tick_spacing = 2
fig = plt.figure(figsize=(8, 6))
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(5,3),
                 axes_pad=0.0,
                 aspect=False
                 )
# Read in user input to set the patch, blind, zmin,zmax, nbootstrap
if len(sys.argv) <2: 
    print("Usage: %s LFVER e.g 2Dbins_5x50" % sys.argv[0])
    sys.exit(1)
else:
	LFVER=sys.argv[1] # catalogue version identifier

# number of tomographic bins, and band power modes to plot
ntomobin=5
ntomocomb=15
nmodes=8

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
#covdat='/disk05/calin/91_Data/mockFootprint/zosterops/MFP_for_others/cov_sim_PeeB_obs.dat'
#These are 3x2pt covs, but the first 120x120 should be the cosmic shear band powers
#Bmode covariance for the null-test
Bcovdat='//home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/Calc_2pt_Stats/PLOTTING/methodology_paper_cov/thps_cov_kids1000_egretta_bp_apo_obs_bandpower_B_apod_matrix.dat'
Bcov=np.loadtxt(Bcovdat)
#Emode covariance to compare the amplitude of the Bmode to the expected Emode signal
Ecovdat='//home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/Calc_2pt_Stats/PLOTTING/methodology_paper_cov/thps_cov_kids1000_egretta_bp_apo_obs_bandpower_E_apod_matrix.dat'
Ecov=np.loadtxt(Ecovdat)

# a guess at the SOM neff difference
#Bcov = Bcov*1.2*1.2

#but the first 8x3 rows are w (2 bins, auto, auto cross)
# the next 8x5x2 rows are gammat (8 nodes, 5 sources, 2 lenses)
# then it should be the 8x15 cosmic shear band powers  
startpt=nmodes*(3+10)
# and set up a smaller array for each tomobin combination
covizjz=np.zeros((nmodes,nmodes))
# and an array for the cosmic shear only components
cov_cs=Bcov[startpt:startpt+nmodes*ntomocomb,startpt:startpt+nmodes*ntomocomb]

# tomake different sub plots we count the grid square that we want to plot in
#initialising the counter
gridpos=-1

#information about the file names
filetop='/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/Pkk/xi2bandpow_output_K1000_ALL_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c'
filetail='nbins_8_Ell_100.0_1500.0_zbins'

# theory curves
#read in the expectation value for the Emode cosmic shear signal
MD='/home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/'

# read in B mode data per tomo bin combination
for iz in range(1,ntomobin+1):
    for jz in range(iz,ntomobin+1):

        gridpos=gridpos + 1
        ax=grid[gridpos]

        # read in the data
        tomochar='%s_%s'%(iz,jz)
        Bnfile='%s_%s_%s_%s.dat'%(filetop,LFVER,filetail,tomochar)
        indata = np.loadtxt(Bnfile)
        ell=indata[:,0]
        PBB=indata[:,3]  # this is l^2 P_BB/ 2pi
        #PBBerr=indata[:,4] # a guess at the error, but better to use the mock or full analytical covariance

        #the expected bandpower_shear_e are l^2 Cl_E/2pi
        BPtheory=np.loadtxt('%s/ForBG/outputs/test_output_S8_fid_test/bandpower_shear_e/bin_%d_%d.txt'%(MD,jz,iz))
        ellmin=np.loadtxt('%s/ForBG/outputs/test_output_S8_fid_test/bandpower_shear_e/l_min_vec.txt'%(MD))
        ellmax=np.loadtxt('%s/ForBG/outputs/test_output_S8_fid_test/bandpower_shear_e/l_max_vec.txt'%(MD))
        elltheory = (ellmax+ellmin)*0.5
        #ax.plot(elltheory,BPtheory/elltheory*1e7,color='blue',label='$\\P_E/100$')

        #breaking up the large Emode covariance matrix to plot the Emode significance for this
        # tomographic bin combination                
        ipos=startpt+gridpos*nmodes
        cov_izjz=Ecov[ipos:ipos+nmodes,ipos:ipos+nmodes]
        diagerr=np.sqrt(np.diagonal(cov_izjz))

        BP_high = (BPtheory + diagerr)/elltheory*1e7
        BP_low = (BPtheory - diagerr)/elltheory*1e7
        ax.fill_between(elltheory, BP_low, BP_high, color='lightgrey',label='$P_E$')
        
        #breaking up the large Bmode covariance matrix to find the significance
        #of the B mode for this tomographic bin combination
        ipos=startpt+gridpos*nmodes
        cov_izjz=Bcov[ipos:ipos+nmodes,ipos:ipos+nmodes]

        if iz==2 :
            print(jz,np.corrcoef(cov_izjz))
        
        diagerr=np.sqrt(np.diagonal(cov_izjz))
        invcov=np.linalg.inv(cov_izjz)
        # calculate the null chi-sq value and associated p-value
        chisq_null = np.matmul(PBB,(np.matmul(invcov,PBB)))
        pval=chi2.sf(chisq_null,nmodes)
        pvalchar='p=%0.3f'%(pval)

        # now plot the results (present l PBB/2pi rather than l^2 PBB/2pi which is given in the data file) 
        # inclue with annotations of the bin combination and p-value 
        ax.errorbar(ell, PBB/ell*1e7, yerr=diagerr/ell*1e7, fmt='o', color='magenta',label=tomochar,markerfacecolor='none')
        ax.axhline(y=0, color='black', ls=':')
        ax.set_ylim(-2.8,5.9)
        ax.annotate(tomochar, xy=(0.2,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.annotate(pvalchar, xy=(0.95,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        ax.set_xscale('log')
        ax.set_xlim(101,1700.0)

        # I also want to know the significance over all bin combinations
        # to do this I need to construct a large data vector

        if gridpos==0:
            PBB_all = PBB
        else:
            PBB_all=np.append(PBB_prev,PBB)
        PBB_prev = PBB_all

# calculate the chi-sq null test for the full data vector
invcov=np.linalg.inv(cov_cs)
chisq_null = np.matmul(PBB_all,(np.matmul(invcov,PBB_all)))
pval=chi2.sf(chisq_null,nmodes*ntomocomb)
# write out the full p-value to the plot title truncated to 3dp
titlechar='%s p=%0.3f'%(LFVER,pval)
grid[1].set_title(titlechar)
#write out to screen the full p-value, rather than the 3dp truncated version
print('p=%5.3e, chi-sq=%8.2f'%(pval,chisq_null))

#add labels
grid[6].set_ylabel('$\ell P_B / 2\pi \,\, [10^{-7}]$')
grid[12].set_xlabel('$\ell$')
grid[13].set_xlabel('$\ell$')
grid[14].set_xlabel('$\ell$')

outfile='figures/PBB_Bmodes_%s.png'%(LFVER)
plt.savefig(outfile)
plt.show()


