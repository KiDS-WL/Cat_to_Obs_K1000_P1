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

# number of tomographic bins, COSEBIS modes and the value of the modes
ntomobin=5
nmodes=5
n=np.arange(1, nmodes+1, dtype=int)

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
covdat='/disk05/calin/91_Data/mockFootprint/zosterops/MFP_for_others/cov_sim_Bn_obs.dat'
cov=np.loadtxt(covdat)
# and set up a smaller array for each tomobin combination
covizjz=np.zeros((nmodes,nmodes))

# tomake different sub plots we count the grid square that we want to plot in
#initialising the counter
gridpos=-1

#information about the file names
filetop='/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/COSEBIS/Bn_COSEBIS_K1000_ALL_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c'
filetail='nbins_300_theta_0.5_300_zbins'

# read in B mode data per tomo bin combination
for iz in range(1,ntomobin+1):
    for jz in range(iz,ntomobin+1):

        gridpos=gridpos + 1
        ax=grid[gridpos]

        tomochar='%s_%s'%(iz,jz)
        Bnfile='%s_%s_%s_%s.asc'%(filetop,LFVER,filetail,tomochar)
        Bn = np.loadtxt(Bnfile)

        #breaking up the large covariance matrix to find the significance per tomographic bin combination
        ipos=gridpos*nmodes
        cov_izjz=cov[ipos:ipos+nmodes,ipos:ipos+nmodes]
        diagerr=np.sqrt(np.diagonal(cov_izjz))
        invcov=np.linalg.inv(cov_izjz)
        # calculate the null chi-sq value and associated p-value
        chisq_null = np.matmul(Bn,(np.matmul(invcov,Bn)))
        pval=chi2.sf(chisq_null,nmodes)
        pvalchar='p=%0.3f'%(pval)

        # and plot the results with annotations of the bin combination and p-value
        ax.errorbar(n, Bn*1e10, yerr=diagerr*1e10, fmt='o', color='magenta',label=tomochar,markerfacecolor='none')
        ax.axhline(y=0, color='black', ls=':')
        ax.set_ylim(-9,12)
        ax.annotate(tomochar, xy=(0.2,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.annotate(pvalchar, xy=(0.95,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        # I also want to know the significance over all bin combinations
        # to do this I need to construct a large data vector

        if gridpos==0:
            Bn_all = Bn
        else:
            Bn_all=np.append(Bn_prev,Bn)
        Bn_prev = Bn_all

# calculate the chi-sq null test for the full data vector
invcov=np.linalg.inv(cov)
chisq_null = np.matmul(Bn_all,(np.matmul(invcov,Bn_all)))
pval=chi2.sf(chisq_null,nmodes*15)
# write out the full p-value to the plot title truncated to 3dp
titlechar='%s p=%0.3f'%(LFVER,pval)
grid[1].set_title(titlechar)
#write out to screen the full p-value, rather than the 3dp truncated version
print('p=%5.3e, chi-sq=%8.2f'%(pval,chisq_null))

#add labels
grid[6].set_ylabel('$B_n \, [10^{-10} rad^2]$')
grid[12].set_xlabel('n')
grid[13].set_xlabel('n')
grid[14].set_xlabel('n')

outfile='figures/COSEBIS_Bmodes_%s.png'%(LFVER)
plt.savefig(outfile)
plt.show()


