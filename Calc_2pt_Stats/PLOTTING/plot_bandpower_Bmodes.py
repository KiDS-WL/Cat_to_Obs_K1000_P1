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
nmodes=8

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
covdat='/disk05/calin/91_Data/mockFootprint/zosterops/MFP_for_others/cov_sim_PeeB_obs.dat'
cov=np.loadtxt(covdat)
# and set up a smaller array for each tomobin combination
covizjz=np.zeros((nmodes,nmodes))

# tomake different sub plots we count the grid square that we want to plot in
#initialising the counter
gridpos=-1

#information about the file names
filetop='/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/Pkk/xi2bandpow_output_K1000_ALL_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c'
filetail='nbins_8_Ell_100.0_1500.0_zbins'

# read in B mode data per tomo bin combination
for iz in range(1,ntomobin+1):
    for jz in range(iz,ntomobin+1):

        gridpos=gridpos + 1
        ax=grid[gridpos]

        tomochar='%s_%s'%(iz,jz)
        Bnfile='%s_%s_%s_%s.dat'%(filetop,LFVER,filetail,tomochar)
        indata = np.loadtxt(Bnfile)
        ell=indata[:,0]
        PBB=indata[:,3]  # this is l^2 P_BB/ 2pi
        #PBBerr=indata[:,4] # a guess at the error, but better to use the mock or full analytical covariance

        #breaking up the large covariance matrix to find the significance per tomographic bin combination
        ipos=gridpos*nmodes
        cov_izjz=cov[ipos:ipos+nmodes,ipos:ipos+nmodes]
        diagerr=np.sqrt(np.diagonal(cov_izjz))
        invcov=np.linalg.inv(cov_izjz)
        # calculate the null chi-sq value and associated p-value
        chisq_null = np.matmul(PBB,(np.matmul(invcov,PBB)))
        pval=chi2.sf(chisq_null,nmodes)
        pvalchar='p=%0.3f'%(pval)

        # now plot the results (present l PBB/2pi rather than l^2 PBB/2pi which is given in the data file) 
        # inclue with annotations of the bin combination and p-value 
        ax.errorbar(ell, PBB/ell*1e7, yerr=diagerr/ell*1e7, fmt='o', color='magenta',label=tomochar)
        ax.axhline(y=0, color='black', ls=':')
        ax.set_ylim(-3,3)
        ax.annotate(tomochar, xy=(0.2,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.annotate(pvalchar, xy=(0.95,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        ax.set_xscale('log')

        # I also want to know the significance over all bin combinations
        # to do this I need to construct a large data vector

        if gridpos==0:
            PBB_all = PBB
        else:
            PBB_all=np.append(PBB_prev,PBB)
        PBB_prev = PBB_all

# calculate the chi-sq null test for the full data vector
invcov=np.linalg.inv(cov)
chisq_null = np.matmul(PBB_all,(np.matmul(invcov,PBB_all)))
pval=chi2.sf(chisq_null,nmodes*15)
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


