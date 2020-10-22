# This file contains the most basic code to read in and plot the
# K1000 mass maps. The code to make the mass maps is in Make_K1000_MassMaps.py
# (same directory on cuillin)
# This is just for plotting and pretty-fying. 
# B. Giblin, 24/07/2020

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys

nside = 1024
ps = 60*(4*np.pi * (180./np.pi)**2 / (12*nside**2))**0.5 # pixel scale in arcmin / pxl
sigma = 10. # smoothing scale in arcmin per pixel
INDIR = 'MassMaps' # input directory of the SNR map
SNR = hp.read_map('%s/KALL_SNR0-100_nside%s_sigma%sarcmin_ZBcut0.1-1.2.fits' %(INDIR,nside,sigma))


# Main code to visualise the mass maps
def Plot_Orthview(figsize, Map, rot, mini, maxi, savename, xlimits, ylimits, ras, decs_lo, decs_up):
    fontsize=14
    plt.figure(figsize=figsize)
    hp.orthview(Map, fig=1, min=mini, max=maxi, cbar=None,
                rot=rot, xsize=2*(ras[-1] - ras[0])*60 / ps, # xsize is v important for determining pixelisation
                half_sky=True,cmap='inferno',badcolor='black',
                title=None) #, format='%.1f', flip='geo' )

    # Set x & y limits
    plt.xlim(xlimits)
    plt.ylim(ylimits)

    # ras & decs_up/decs_low are the tick labels to be plotted with:
    for r in ras: # plot the ra labels
        hp.visufunc.projtext(r,decs_up-5,'%.0f'%r+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)                    
    # (might have to play around with ra positioning of dec labels)
    hp.visufunc.projtext(ras[0]-10,decs_lo,'%.0f'%decs_lo+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    hp.visufunc.projtext(ras[0]-10,decs_lo+10,'%.0f'%(decs_lo+10)+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    hp.visufunc.projtext(ras[0]-10,decs_up,'%.0f'%decs_up+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    #hp.visufunc.projplot(ras,decs_up,'wo',lonlat=True) # Plot dots

    # Plots ra/dec contours every dpar intervals in ra, dmer intervals in dec
    hp.graticule(dpar=10, dmer=10, coord='G', color='white')
    plt.savefig(savename)
    plt.show()
    return


# PRESS RELEASE PLOTS
SNR_lo = -5
SNR_hi = 5  # SNR colour map limits

# NORTH
figsize_N = (15,2)
rot_N = [182.,0.]               # rotation applied to centre North patch
xlimits_N = [-0.84,0.82]        # unit value xlimits (range [-1,1])
ylimits_N = [-0.1,0.1]          # unit value y limits (range [-1,1])
ras_N = np.arange(140,220,10)   # where RA labels will appear if 'hp.visufunc.projtext' lines are uncommented
dec_lo_N = -10
dec_hi_N = 10                   # the upper and lower dec limits where dec labels will appear
savename_N = '%s/orthview_KN_SNR0-100_nside%s_sigma%sarcmin_ZBcut0.1-1.2.png' %(INDIR,nside,sigma)

#Plot_Orthview(figsize_N, SNR, rot_N, SNR_lo,SNR_hi,
#              savename_N,
#              xlimits_N, ylimits_N, ras_N, dec_lo_N, dec_hi_N)

print("I've exited the code here, because you have to close Figure 1, containing the N patch map",
      " before making the South patch map. Otherwise they overplot.",
      " There's probably a work-around this but I haven't found it.")
#sys.exit()

# SOUTH
figsize_S = (12,4)              # very roughly in proportion to figsize_S
rot_S = [10.,-30.]
xlimits_S = [-0.65,0.6]        
ylimits_S = [-0.22,0.08]         
ras_S = np.arange(-20,70,10)   
dec_lo_S = -40
dec_hi_S = -20                   
savename_S = '%s/orthview_KS_SNR0-100_nside%s_sigma%sarcmin_ZBcut0.1-1.2.png' %(INDIR,nside,sigma)

# NOTE: use x/y-limits when plotting RA,Dec labels: [-0.65,0.65], [-0.3,0.25]
Plot_Orthview(figsize_S, SNR, rot_S, SNR_lo,SNR_hi,
	      savename_S,
              xlimits_S, ylimits_S, ras_S, dec_lo_S, dec_hi_S)



