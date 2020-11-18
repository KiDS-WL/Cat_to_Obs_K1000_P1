import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from astropy.io import ascii
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

# Some font setting
text_width = 523.5307/72
column_width = 256.0748/72

matplotlib.rc("text", usetex=True)
matplotlib.rc("text.latex", preamble=
r"""
\usepackage{txfonts}
\newcommand{\mathdefault}[1][]{}""")

matplotlib.rc("font", family="Times", size=10)

#set up the figure grid PEE panels, a gap and then the PBB panels

gridspec = dict(hspace=0.0, wspace=0.0, width_ratios=[1, 1])
fig, axes = plt.subplots(nrows=1, ncols=2, gridspec_kw=gridspec, figsize=(text_width, text_width*0.4))

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
#These are 3x2pt covs, even though they are stored in the Pkk_cov directory
#Bcovdat='../Pkk/Pkk_cov/thps_cov_kids1000_mar30_bandpower_B_apod_0_matrix.dat'
#Bcov=np.loadtxt(Bcovdat)

#start the tomographic bin counter
binid=-1

#information about the file names
filetop='wedges_data/BOSS.DR12.'
filetail='z.3xiwedges_measurements.txt'

# theory curves
#read in the expectation value for the Emode cosmic shear signal
#MD='/home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/'
#MD='/Users/macleod/CH_work/Cat_to_Obs_K1000_P1/'
#MD='/Users/heymans/KiDS/Cat_to_Obs_K1000_P1/'
MD='/Users/yooken/Research/KiDS/Cat_to_Obs_K1000_P1'

#Set the x/y limits
xmin=20.0
xmax=160
ymin=-80
ymax=150.0

#which sub plot do we want to plot this in
lens_count=-1

# read in wedges data per lens bin
for lensbin in ("low","high"):

    #increment the subpanel counter
    lens_count = lens_count+1

    # read in the data
    if lensbin=="low":
        labelchar='$0.2 < z \leq 0.5$'
    else:
        labelchar='$0.5 < z \leq 0.75$'

    wedgesfile='%s%s%s'%(filetop,lensbin,filetail)
    indata = np.loadtxt(wedgesfile)
    s=indata[:,0]
    xi1=indata[:,1]  # this is xi transverse
    sigxi1=indata[:,2]  
    xi2=indata[:,3]  # this is xi intermediate
    sigxi2=indata[:,4]  
    xi3=indata[:,5]  # this is xi parallel
    sigxi3=indata[:,6]  
    
    # and read in the expected clustering measurements
    theory=np.loadtxt('%s/Predictions/iterated_cov_MAP_BlindC/xi_wedges/bin_%d.txt'%(MD,lens_count+1))
    stheory = s[4:32] #guess - need to check with Marika that this is the s-scale

    #PLOT THE wedges!
    #which grid cell do we want to plot this in?
    ax=axes[lens_count]
    # only label the subplots at the edges
    ax.label_outer()

    offset=0.5
    #transverse
    trans,=ax.plot(stheory-offset,theory[0]*stheory*stheory,color='m',linewidth=1.5,label="Transverse")
    ax.errorbar(s-offset, xi1*s*s, yerr=sigxi1*s*s, fmt='D', markersize=2.5,color='m')
    #intermediate
    intermediate,=ax.plot(stheory,theory[1]*stheory*stheory,color='blue',linewidth=1.5,label="Intermediate")
    ax.errorbar(s, xi2*s*s, yerr=sigxi1*s*s, fmt='s', markersize=2.5,color='blue',alpha=0.6)
    #parallel
    parallel,=ax.plot(stheory+offset,theory[2]*stheory*stheory,color='black',linewidth=1.5,label="Parallel")
    ax.errorbar(s+offset, xi3*s*s, yerr=sigxi1*s*s, fmt='o', markersize=2.5,color='black',alpha=0.6)

    # set the limits of the plot
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # add the tomographic bin combination
    ax.annotate(labelchar, xy=(0.97,0.95),xycoords='axes fraction',
            size=12, ha='right', va='top')

    # Split legend across both the panels
    if lensbin=="low":
        plt.setp(ax.get_xticklabels()[-1], visible=False)
        ax.legend([trans,intermediate],["Transverse","Intermediate"],loc='lower left',fontsize=10, frameon=False)
    else:
        ax.legend([parallel],["Parallel"],loc='lower left',fontsize=10, frameon=False)

#add plot labels
axes[0].set_xlabel('s  [$h^{-1}$ Mpc]')
axes[1].set_xlabel('s  [$h^{-1}$ Mpc]')
axes[0].set_ylabel('s$^2 \, \\xi_{gg}$ [$h^{-2}$ Mpc$^2$]')

plt.tight_layout()

outfile='BOSS_Sanchez_wedges.pdf'
plt.savefig(outfile,dpi=300)
#plt.show()


