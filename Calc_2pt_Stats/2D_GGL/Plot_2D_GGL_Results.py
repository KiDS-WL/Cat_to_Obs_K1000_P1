# 04/09/2019, B.M. Giblin, PhD Student (sort of), Edinburgh
# Plot the 2D gamma_t, gamma_x results from Run_GGL_K1000.py


import numpy as np
import treecorr
import sys
from astropy.io import fits
import os
import time

# Plotting modules
from matplotlib import rcParams
import pylab as plt
import matplotlib.gridspec as gridspec

#rc('text',usetex=True) #cuillin is v unhappy with this command?
#rc('font',size=14)
#rc('legend',**{'fontsize':14})
#rc('font',**{'family':'serif','serif':['Computer Modern']})

rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}
plt.rc('font', **font)

#This makes episilons appear as epsilons rather than varepsilons
plt.rcParams["mathtext.fontset"] = "cm"


NorS = sys.argv[1]  # N for North, S for South
Blind = sys.argv[2] # Just put 'A' for now

ZBlo = sys.argv[3]  # If numbers are inputted, this Z_B cut is applied to sources
ZBhi = sys.argv[4]  # For no ZBcut, put None for both of these

Zlo = sys.argv[5]  # If numbers are inputted, this Z cut is applied to lenses and randoms
Zhi = sys.argv[6]  # For no ZBcut, put None for both of these

# Sources
if str(ZBlo) == "None":
    ZBlabel = "None"
else:
    ZBlabel = "%s-%s" %(ZBlo, ZBhi)

# Lenses/randomsls 
if str(Zlo) == "None":
    Zlabel = "None"
else:
    Zlabel = "%s-%s" %(Zlo, Zhi)

# Read in the 2D gamma_t, gamma_x results
Ntiles = 475
dxdy = np.load('Results/K%s.Blind%s.2Ddxdy.Ntiles%s.SourceZBcut%s.LensZcut%s.npy'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))
# with lenses
gamma_t_2D = np.load('Results/K%s.Blind%s.2Dgamma_t.Ntiles%s.SourceZBcut%s.LensZcut%s.npy'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))
gamma_x_2D = np.load('Results/K%s.Blind%s.2Dgamma_x.Ntiles%s.SourceZBcut%s.LensZcut%s.npy'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))
# with randoms
gamma_tr_2D = np.load('Results/K%s.Blind%s.2Dgamma_t_randoms.Ntiles%s.SourceZBcut%s.LensZcut%s.npy'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))
gamma_xr_2D = np.load('Results/K%s.Blind%s.2Dgamma_x_randoms.Ntiles%s.SourceZBcut%s.LensZcut%s.npy'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))


# Read the 1D gamma_t binned in theta 
# with lenses
theta, gamma_t_theta, gamma_x_theta = np.loadtxt('Results/K%s.Blind%s.gamma_tx.SourceZBcut%s.LensZcut%s.dat'%(NorS,Blind,ZBlabel,Zlabel), usecols=(0,2,3), unpack=True)
# with randoms
theta, gamma_tr_theta, gamma_xr_theta = np.loadtxt('Results/K%s.Blind%s.gamma_tx_randoms.SourceZBcut%s.LensZcut%s.dat'%(NorS,Blind,ZBlabel,Zlabel), usecols=(0,2,3), unpack=True)

# Sort gamma_t_theta onto a symmetrical grid 
gamma_t_1D = np.zeros_like( gamma_t_2D )
gamma_x_1D = np.zeros_like( gamma_t_2D )
gamma_tr_1D = np.zeros_like( gamma_t_2D )
gamma_xr_1D = np.zeros_like( gamma_t_2D )
ps = 0.213 	# arcsec per pxl

dxdy_arcmin = dxdy * ps / 60
nbins = gamma_t_2D.shape[0]
for i in range(nbins):
	for j in range(nbins):
		gamma_t_1D[i,j] = np.interp( dxdy_arcmin[i,j], theta, gamma_t_theta )
		gamma_tr_1D[i,j] = np.interp( dxdy_arcmin[i,j], theta, gamma_tr_theta )
		gamma_x_1D[i,j] = np.interp( dxdy_arcmin[i,j], theta, gamma_x_theta )
		gamma_xr_1D[i,j] = np.interp( dxdy_arcmin[i,j], theta, gamma_xr_theta )


# Get the angular values for plotting on the x/y axes
# NOTE: THIS WILL ONLY BE EXACT IF ARRAY HAS AN ODD NUMBER OF BINS
if dxdy.shape[0] % 2 == 0:
	print("WARNING: Your 2D array has an even number of bins: %sx%s. This means angular values plotted on the x,y axes of following plot will not be quite right." %(dxdy.shape[0],dxdy.shape[1]) )
	print("It's best to recompute the array with odd number of bins on each side.")
nrow = np.int(dxdy.shape[0] / 2)	# This will give the row on level with the origin
                                        # (i.e. row 10 for a 21x21 sized array)
# Make the angular x/y-values array, symmetric about the origin.
# Checked that you get identical results if you calculate them from cols instead of rows.
ang_vals = np.copy(dxdy_arcmin[nrow,:])
ang_vals[ang_vals<1e-10] = 0
ang_vals[:nrow] *= -1


# Plot the gamma_{t,x} 2D and 1D
cmap = plt.cm.rainbow
f = plt.figure(figsize = (10,24))
gs1 = gridspec.GridSpec(4, 2)

#vmin = -0.002 	#gamma_t_2D.min()
#vmax = 0.01    #gamma_t_2D.max()
vmin = -0.0075 	#gamma_t_2D.min()
vmax = 0.035    #gamma_t_2D.max()


# ----------- ROW 1 -----------
# gamma_t^2D w/ lenses
ax1 = plt.subplot(gs1[0], adjustable='box')
cplot = ax1.imshow(gamma_t_2D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]])

ax1.set_ylabel(r'$\Delta\theta_y \, [\rm{arcmin}]$')
ax1.set_title(r'$\epsilon_{\rm t}^{\rm{obs}}$')

# gamma_x^2D w/ lenses
ax2 = plt.subplot(gs1[1], adjustable='box')
cplot = ax2.imshow(gamma_x_2D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]])  
ax2.set_ylabel(r'BOSS$_{2D}$'+'\n')
ax2.yaxis.tick_right()
ax2.set_title(r'$\epsilon_\times^{\rm{obs}}$')


# ----------- ROW 3 -----------
# gamma_t^1D w/ lenses
ax3 = plt.subplot(gs1[2], adjustable='box')
cplot = ax3.imshow(gamma_t_1D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]])
ax3.set_ylabel(r'$\Delta\theta_y \, [\rm{arcmin}]$')

# gamma_x^1D w/ lenses
ax4 = plt.subplot(gs1[3], adjustable='box')
cplot = ax4.imshow(gamma_x_1D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]]) 
ax4.set_ylabel(r'BOSS$_{1D}$'+'\n')
ax4.yaxis.tick_right()

# ----------- ROW 3 -----------
# Delta gamma_t
ax5 = plt.subplot(gs1[4], adjustable='box')
#Delta_gamma_t = gamma_t_2D - gamma_tr_2D - gamma_t_1D
Delta_gamma_t = gamma_t_2D - gamma_t_1D
cplot = ax5.imshow(Delta_gamma_t, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]]) 
ax5.set_ylabel(r'$\Delta\theta_y \, [\rm{arcmin}]$')

# Delta gamma_x
ax6 = plt.subplot(gs1[5], adjustable='box')
#Delta_gamma_x = gamma_x_2D - gamma_xr_2D - gamma_x_1D
Delta_gamma_x = gamma_x_2D - gamma_x_1D
cplot = ax6.imshow(Delta_gamma_x, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]]) 
ax6.set_ylabel(r'BOSS$_{2D - 1D}$' + '\n')
ax6.yaxis.tick_right()

# ----------- ROW 4 -----------
# gamma_t^2D w/ randoms
ax7 = plt.subplot(gs1[6], adjustable='box')
cplot = ax7.imshow(gamma_tr_2D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]]) 
ax7.set_ylabel(r'$\Delta\theta_y \, [\rm{arcmin}]$')
ax7.set_xlabel(r'$\Delta\theta_x \, [\rm{arcmin}]$')

# gamma_x^2D w/ randoms
ax8 = plt.subplot(gs1[7], adjustable='box')
cplot = ax8.imshow(gamma_xr_2D, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, 
					origin='lower', extent=[ang_vals[0],ang_vals[-1],ang_vals[0],ang_vals[-1]]) 
ax8.set_ylabel(r'randoms$_{2D}$'+'\n')
ax8.yaxis.tick_right()
ax8.set_xlabel(r'$\Delta\theta_x \, [\rm{arcmin}]$')


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax6.get_xticklabels(), visible=False)

plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)

f.subplots_adjust(hspace=0, wspace=0)
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.8, 0.15, 0.03, 0.7])
clb=f.colorbar(cplot, cax=cbar_ax) 
clb.ax.set_title('$\epsilon$')

plt.savefig('Results/Plot_K%s.Blind%s.2Dgamma_tx.Ntiles%s.SourceZBcut%s.LensZcut%s.png'%(NorS,Blind,Ntiles,ZBlabel,Zlabel))
plt.show()







