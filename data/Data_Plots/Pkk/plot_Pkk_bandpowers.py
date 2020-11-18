import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import matplotlib
import sys
from astropy.io import ascii
import matplotlib.ticker as ticker

#TT's pre-amble
text_width = 523.5307/72
column_width = 256.0748/72

matplotlib.rc("text", usetex=True)
matplotlib.rc("text.latex", preamble=
r"""
\usepackage{txfonts}
\newcommand{\mathdefault}[1][]{}""")

matplotlib.rc("font", family="Times", size=10)

#set up the figure grid
tick_spacing = 2
fig,axes= plt.subplots(6,7,figsize=(text_width, text_width*0.75),gridspec_kw={'hspace': 0, 'wspace': 0})

# Read in user input to set the patch, blind, zmin,zmax, nbootstrap
#if len(sys.argv) <2: 
#    print("Usage: %s LFVER BLIND e.g 2Dbins_v2_goldclasses_Flag_SOM_Fid A" % sys.argv[0])
#    sys.exit(1)
#else:
#    LFVER=sys.argv[1] # catalogue version identifier
#    BLIND=sys.argv[2] # blind

LFVER="2Dbins_v2_goldclasses_Flag_SOM_Fid" 
BLIND="C"

# number of tomographic bins, and band power modes to plot
ntomobin=5
ntomocomb=15
nmodes=8

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
#These are 3x2pt covs
Bcovdat='Pkk_cov/thps_cov_kids1000_mar30_bandpower_B_apod_0_matrix.dat'
Bcov=np.loadtxt(Bcovdat)
#Emode covariance to compare the amplitude of the Bmode to the expected Emode signal
Ecovdat='Pkk_cov/thps_cov_kids1000_mar30_bandpower_E_apod_0_matrix.dat'
Ecov=np.loadtxt(Ecovdat)


# The first 8x2 rows are w (11, 12)
# The next 8x5 rows are GGL (11, 12, 13, 14, 15) 
# The next 8 rows are w (22)
# Then 8x5 rows for GGL again (21, 22, 23, 24, 25)
# then it there are the 8x15 cosmic shear band powers
startpt=nmodes*(3+10)
# and set up a smaller array for each tomobin combination
covizjz=np.zeros((nmodes,nmodes))

#start the tomographic bin counter
binid=-1

#information about the file names
filetop='Pkk_data/xi2bandpow_output_K1000_ALL_BLIND_'+str(BLIND)+'_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c'
filetail='nbins_8_Ell_100.0_1500.0_zbins'

# theory curves
#read in the expectation value for the Emode cosmic shear signal
#MD='/home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/'
#MD='/Users/macleod/CH_work/Cat_to_Obs_K1000_P1/'
#MD='/Users/heymans/KiDS/Cat_to_Obs_K1000_P1/'
MD='/Users/yooken/Research/KiDS/Cat_to_Obs_K1000_P1'

#Set the x/y limits
xmin=101
xmax=1700.0
ymin=-4.0
Eymax=10.1
Bymax=5.1

# read in B mode data per tomo bin combination
for iz in range(1,ntomobin+1):
    for jz in range(iz,ntomobin+1):

        #increment the bin counter
        binid=binid+1

        # read in the data
        tomochar='%s-%s'%(iz,jz)
        EBnfile='%s_%s_%s_%s.dat'%(filetop,LFVER,filetail,tomochar.replace("-","_"))
        indata = np.loadtxt(EBnfile)
        ell=indata[:,0]
        PBB=indata[:,3]  # this is l^2 P_BB/ 2pi
        PEE=indata[:,1]  # this is l^2 P_EE/ 2pi
        
        # breaking up the large E/Bmode covariance matrix to plot the E/Bmode significance for this
        # tomographic bin combination                
        ipos=startpt+binid*nmodes
        # E mode
        Ecov_izjz=Ecov[ipos:ipos+nmodes,ipos:ipos+nmodes]
        Ediagerr=np.sqrt(np.diagonal(Ecov_izjz))
        # B mode
        Bcov_izjz=Bcov[ipos:ipos+nmodes,ipos:ipos+nmodes]
        Bdiagerr=np.sqrt(np.diagonal(Bcov_izjz))
        
        BPtheory=np.loadtxt('%s/Predictions/iterated_cov_MAP_BlindC/bandpower_shear_e/bin_%d_%d.txt'%(MD,jz,iz))
        ellmin=np.loadtxt('%s/Predictions/iterated_cov_MAP_BlindC/bandpower_shear_e/l_min_vec.txt'%(MD))
        ellmax=np.loadtxt('%s/Predictions/iterated_cov_MAP_BlindC/bandpower_shear_e/l_max_vec.txt'%(MD))
        elltheory = ell #(ellmax+ellmin)*0.5

        #PLOT THE EMODES!
        #which grid cell do we want to plot this in?
        grid_x_E=(iz-1)
        grid_y_E=(5-jz)
        ax=axes[grid_y_E,grid_x_E]

        ax.plot(elltheory,BPtheory/elltheory*1e7,color='m',linewidth=1.5)
        ax.errorbar(ell, PEE/ell*1e7, yerr=Ediagerr/ell*1e7, fmt='o', markersize=2.5,color='b')

        # only label the subplots at the edges
        ax.label_outer()
        # and then override for the autocorrelations
        if iz==jz:
            ax.xaxis.tick_bottom()
            if iz==1:    
                ax.set_xlabel('$\ell$')

        # set the limits of the plot
        ax.set_xscale('log')
        ax.set_xlim(xmin,xmax)
        #ax.set_yscale('log')
        #ax.set_ylim(0.1,160.0)
        ax.set_ylim(ymin,Eymax)

        # add the tomographic bin combination
        ax.annotate(tomochar, xy=(0.95,0.9),xycoords='axes fraction',
            size=10, ha='right', va='top')
        ax.axhline(y=0, color='black', ls=':')
        
        #PLOT THE BMODES!

        #which grid cell do we want to plot this in?
        grid_x_B=(jz+1)
        grid_y_B=(6-iz)
        ax=axes[grid_y_B,grid_x_B]
        ax.errorbar(ell, PBB/ell*1e7, yerr=Bdiagerr/ell*1e7, fmt='o', markersize=2.5,color='b')
        
        # only label the subplots at the edges
        ax.label_outer()
        # and then override for the bin 5 left column
        if jz==5:
            ax.yaxis.tick_right()
        if ((iz==1) and (jz==1)):    
            ax.set_xlabel('$\ell$')

        # set the limits of the plot
        ax.set_ylim(ymin,Bymax)
        ax.set_xscale('log')
        ax.set_xlim(xmin,xmax)

        # add the tomographic bin combination and a horizontal line
        ax.annotate(tomochar, xy=(0.95,0.9),xycoords='axes fraction',
                    size=10, ha='right', va='top')
        ax.axhline(y=0, color='black', ls=':')


#add plot labels
axes[2,0].set_ylabel('$\ell C_{EE} \, / \,  2\pi \,\, [10^{-7}]$')
axes[3,6].yaxis.label.set_visible(True)
axes[3,6].set_ylabel('$\ell C_{BB} \, /\,  2\pi \,\, [10^{-7}]$')
axes[3,6].yaxis.set_label_position('right')

#Finally Blank out the empty cells
for i in range(6):
    blankgrid=6-i
    axes[i,blankgrid].set_visible(False)
    axes[i,blankgrid-1].set_visible(False)

plt.tight_layout()

outfile='Pkk_K1000_%s_%s.pdf'%(LFVER,BLIND)
plt.savefig(outfile,dpi=300)
plt.show()


