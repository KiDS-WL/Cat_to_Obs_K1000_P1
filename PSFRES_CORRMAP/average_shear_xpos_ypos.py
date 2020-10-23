import numpy as np
import ldac
import matplotlib.pyplot as plt

import sys
from scipy.stats import binned_statistic_2d

from matplotlib import rcParams
from astropy.io import fits
from fitting import * # MV scripts
from scipy.stats import binned_statistic

plt.rcParams.update({'font.size': 15})

# standard error of the mean
def SEM(x):
        return np.std(x)/len(x)

# straight line y = f(x)
def func(x,a,b): 
        return a*x+b

md=sys.argv[1]
patch=sys.argv[2]
infile=sys.argv[3]

# Define custom colormaps: Set pixels with no sources to white
cmap = plt.cm.seismic
cmap.set_bad('w', 1.)

##########################################################
# read in the c1 and c2 maps
hdu_list_c1 = fits.open(md+'/c1_map.fits')
c1_map_tmp = hdu_list_c1[0].data

hdu_list_c2 = fits.open(md+'/c2_map.fits')
c2_map_tmp = hdu_list_c2[0].data

# How big is this array so we can bin the data with the same dimensions?
Npix_fine=len(c1_map_tmp)
# Rebinned size of array dimension:
Npix=70

# Rebin
c1_map_rebin = c1_map_tmp.reshape(Npix, Npix_fine/Npix, Npix, Npix_fine/Npix)
c1_map = c1_map_rebin.mean(axis=3).mean(axis=1)

c2_map_rebin = c2_map_tmp.reshape(Npix, Npix_fine/Npix, Npix, Npix_fine/Npix)
c2_map = c2_map_rebin.mean(axis=3).mean(axis=1)

# find the edges of the image from the exposure map
hdu_list_exp = fits.open(md+'/exposure_map.fits')
exp_map_tmp = hdu_list_exp[0].data
exp_map_rebin = exp_map_tmp.reshape(Npix, Npix_fine/Npix, Npix, Npix_fine/Npix)
exp_map = exp_map_rebin.mean(axis=3).mean(axis=1)

# and clip so it is zero for out of the image and 1 within
w_map = np.clip(exp_map,0,1)
w_map = w_map.astype(int)
#w_map[61:,:] = 0
#w_map[:10,:] = 0
#w_map[:,61:] = 0
#w_map[:,:10] = 0

max = np.max(w_map)
min = np.min(w_map)
max = np.max((np.abs(max), np.abs(min)))
plt.imshow(w_map, origin='lower',
           aspect='equal', interpolation='none',
           cmap=cmap, clim=(-max, max)) #, vmin=-0.005, vmax=0.005)
plt.colorbar()
plt.savefig(md+'/PLOTS/w_map_%s_V1.0.0A.png' % (patch))
plt.close()

# remove average c1/c2
dec1 = np.sum(w_map*c1_map)/np.sum(w_map)
dec2 = np.sum(w_map*c2_map)/np.sum(w_map)

c1_map = (c1_map - dec1)*w_map
c2_map = (c2_map - dec2)*w_map

max = np.max(c1_map)
min = np.min(c1_map)
max = np.max((np.abs(max), np.abs(min)))
plt.imshow(c1_map, origin='lower',
           aspect='equal', interpolation='none',
           cmap=cmap, clim=(-max, max), extent=[0,21000,0,21000]) #, vmin=-0.005, vmax=0.005)
plt.xlabel("$x$ [pix]", labelpad=-1)
plt.ylabel("$y$ [pix]") #, labelpad=-2)
cbar = plt.colorbar()
cbar.set_label('$c_1$',size=18)
plt.savefig(md+'/PLOTS/c1_map_%s_V1.0.0A.pdf' % (patch))
plt.close()

max = np.max(c2_map)
min = np.min(c2_map)
max = np.max((np.abs(max), np.abs(min)))
plt.imshow(c2_map, origin='lower',
           aspect='equal', interpolation='none',
           cmap=cmap, clim=(-max, max)) #, vmin=-0.005, vmax=0.005)
plt.colorbar()
plt.savefig(md+'/PLOTS/c2_map_%s_V1.0.0A.pdf' % (patch))
plt.close()

###########################################################
#open the ldac catalogue using functions in ldac.py
ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the useful columns
KiDS_RA=ldac_table['ALPHA_J2000']
KiDS_Dec=ldac_table['DELTA_J2000']
Xpos_in=ldac_table['Xpos']
Ypos_in=ldac_table['Ypos']

e1_in=ldac_table['autocal_e1_A']
e2_in=ldac_table['autocal_e2_A']
weight_in=ldac_table['recal_weight_A']
mag_in=ldac_table['MAG_AUTO']
mask=ldac_table['MASK']

we1_in=e1_in*weight_in
we2_in=e2_in*weight_in

# make a magnitude selection
magsel = (mag_in>20.)

# make MASK selection
masksel = (mask==0)

# combine selections
allsel = magsel * masksel

Xpos = Xpos_in[allsel]
Ypos = Ypos_in[allsel]
we1 = we1_in[allsel]
we2 = we2_in[allsel]
weight = weight_in[allsel]
mag = mag_in[allsel]

c1 = np.sum(we1)/np.sum(weight)
c2 = np.sum(we2)/np.sum(weight)
print c1, c2
##### Deliberately introduce a de term ###
#Xpos_int = (Xpos/(21000/Npix)).astype(int)
#Ypos_int = (Ypos/(21000/Npix)).astype(int)
#beta1_fake = 0. #10. 
#beta2_fake = 0. #-20. 
#we1 = we1 + beta1_fake*c1_map[Ypos_int,Xpos_int]*weight
#we2 = we2 + beta2_fake*c2_map[Ypos_int,Xpos_int]*weight
##########################################

# we want the binning of the data to be the same as the binning of the c_maps so
ibn=Npix

# binned_statistic has it's axes flipped - trial and error found this way of ordering
# puts the resulting array on the same axes as the cmap

e1_sum, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, we1,    statistic='sum', bins=ibn, range=[[0.0,21000.0],[0.0,21000.0]])
e2_sum, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, we2,    statistic='sum', bins=ibn, range=[[0.0,21000.0],[0.0,21000.0]])
w_sum,  yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, weight, statistic='sum', bins=ibn, range=[[0.0,21000.0],[0.0,21000.0]])
w_sum_filter = w_sum>0

w_sum_flat = w_sum.flatten()
w_sum_flat_filter = w_sum_flat>0

# reject pixels with zero weight
# and subtract off the average c-term only within the image
e1_mean = np.zeros((Npix,Npix))
e2_mean = np.zeros((Npix,Npix))
e1_mean[w_sum_filter] = (e1_sum[w_sum_filter]/w_sum[w_sum_filter])*w_map[w_sum_filter] - c1*w_map[w_sum_filter]
e2_mean[w_sum_filter] = (e2_sum[w_sum_filter]/w_sum[w_sum_filter])*w_map[w_sum_filter] - c2*w_map[w_sum_filter]

max = np.max(e1_mean)
min = np.min(e1_mean)
max = np.max((np.abs(max), np.abs(min))) / 2.
plt.imshow(e1_mean, origin='lower',
           aspect='equal', interpolation='none',
           cmap=cmap, clim=(-max, max)) #, vmin=-0.005, vmax=0.005)
plt.colorbar()
plt.savefig(md+'/PLOTS/e1_mean_%s_V1.0.0A.png' % (patch))
plt.close()

max = np.max(e2_mean)
min = np.min(e2_mean)
max = np.max((np.abs(max), np.abs(min))) / 2.
plt.imshow(e2_mean, origin='lower',
           aspect='equal', interpolation='none',
           cmap=cmap, clim=(-max, max)) #, vmin=-0.005, vmax=0.005)
plt.colorbar()
plt.savefig(md+'/PLOTS/e2_mean_%s_V1.0.0A.png' % (patch))
plt.close()

### bin the data points ###
minx = np.min(c1_map.flatten())
maxx = np.max(c1_map.flatten())
stepsize1 = (maxx-minx)/10.
bins = np.arange(start=minx, stop=maxx+stepsize1, step=stepsize1)
bin_centres1 = (bins[:-1]+bins[1:])/2
c1_map_binned  = binned_statistic(c1_map.flatten()[w_sum_flat_filter], c1_map.flatten()[w_sum_flat_filter], bins=bins)[0]
e1_mean_binned = binned_statistic(c1_map.flatten()[w_sum_flat_filter], e1_mean.flatten()[w_sum_flat_filter], bins=bins)[0]
e1_SEM_binned  = binned_statistic(c1_map.flatten()[w_sum_flat_filter], e1_mean.flatten()[w_sum_flat_filter], statistic=SEM, bins=bins)[0]

minx = np.min(c2_map.flatten())
maxx = np.max(c2_map.flatten())
stepsize2 = (maxx-minx)/10.
bins = np.arange(start=minx, stop=maxx+stepsize2, step=stepsize2)
bin_centres2 = (bins[:-1]+bins[1:])/2
c2_map_binned  = binned_statistic(c2_map.flatten()[w_sum_flat_filter], c2_map.flatten()[w_sum_flat_filter], bins=bins)[0]
e2_mean_binned = binned_statistic(c2_map.flatten()[w_sum_flat_filter], e2_mean.flatten()[w_sum_flat_filter], bins=bins)[0]
e2_SEM_binned  = binned_statistic(c2_map.flatten()[w_sum_flat_filter], e2_mean.flatten()[w_sum_flat_filter], statistic=SEM, bins=bins)[0]
###########################

# now we need to find the best fit of the c_map to the e_map
# flattern the 2D arrays into 1D for the fit
# calculate alpha and c for the full sample
p0=[0.,0.]
bounds=((-100.,-100.),(100.,100.))

param, err, rchi, dof, cov =general_fit(
        f=func,
        xdata=c1_map.flatten()[w_sum_flat_filter],
        ydata=e1_mean.flatten()[w_sum_flat_filter],
        p0=p0,
        sigma=1.0/w_sum_flat[w_sum_flat_filter],
        bounds=bounds,
        method='trf',
        loss='soft_l1',
        f_scale=0.5,
        verbose=1,
        jac='3-point',
        max_nfev=100000)
print "Fit for c1"
print param, err, rchi, dof
print cov
print

#plt.errorbar(c1_map.flatten()[w_sum_flat_filter], e1_mean.flatten()[w_sum_flat_filter], fmt='x', yerr=1.0/w_sum_flat[w_sum_flat_filter], alpha=0.5)
plt.plot(c1_map.flatten()[w_sum_flat_filter], e1_mean.flatten()[w_sum_flat_filter], ',', alpha=0.2)
#plt.errorbar(c1_map_binned, e1_mean_binned, yerr=e1_SEM_binned, fmt='o')
plt.errorbar(bin_centres1, e1_mean_binned, yerr=e1_SEM_binned, xerr=0.5*stepsize1, fmt='o')
minx=np.min(c1_map.flatten())
maxx=np.max(c1_map.flatten())
plt.plot((minx,maxx),(param[0]*minx+param[1],param[0]*maxx+param[1]))
#plt.plot((minx,maxx),(beta1_fake*minx,param[0]*maxx))
plt.xlabel('c1_map', labelpad=-1)
plt.ylabel('e1_mean')
plt.ylim(-0.01,0.01)
plt.savefig(md+'/PLOTS/c1_e1_%s_V1.0.0A.png' % (patch))
plt.close()

param, err, rchi, dof, cov =general_fit(
        f=func,
        xdata=c2_map.flatten()[w_sum_flat_filter],
        ydata=e2_mean.flatten()[w_sum_flat_filter],
        p0=p0,
        sigma=1.0/w_sum_flat[w_sum_flat_filter],
        bounds=bounds,
        method='trf',
        loss='soft_l1',
        f_scale=0.5,
        verbose=1,
        jac='3-point',
        max_nfev=100000)
print "Fit for c2"
print param, err, rchi, dof
print cov
print

#plt.errorbar(c2_map.flatten()[w_sum_flat_filter], e2_mean.flatten()[w_sum_flat_filter], fmt='x', yerr=1.0/w_sum_flat[w_sum_flat_filter], alpha=0.5)
plt.plot(c2_map.flatten()[w_sum_flat_filter], e2_mean.flatten()[w_sum_flat_filter], ',', alpha=0.2)
#plt.errorbar(c2_map_binned, e2_mean_binned, yerr=e2_SEM_binned, fmt='o')
plt.errorbar(bin_centres2, e2_mean_binned, yerr=e2_SEM_binned, xerr=0.5*stepsize2, fmt='o')
minx=np.min(c2_map.flatten())
maxx=np.max(c2_map.flatten())
plt.plot((minx,maxx),(param[0]*minx+param[1],param[0]*maxx+param[1]))
#plt.plot((minx,maxx),(beta2_fake*minx,param[0]*maxx))
plt.xlabel('c2_map', labelpad=-1)
plt.ylabel('e2_mean')
plt.ylim(-0.01,0.01)
plt.savefig(md+'/PLOTS/c2_e2_%s_V1.0.0A.png' % (patch))
plt.clf()

####for i in range(1,2):
####    
####    if(i==1):
####        plt.subplot(121)
####        bindat = e1_mean
####    else:
####        plt.subplot(122)
####        bindat = e2_mean
####
####        #im = plt.imshow(bindat, origin='lower',
####        #                extent=[xedges[-1], xedges[0],
####        #                        yedges[-1], yedges[0]],
####        #                aspect='auto', interpolation='none', cmap=cmap)
####    im = plt.imshow(bindat, origin='lower',
####                    extent=[xedges[0], xedges[ibn],
####                            yedges[0], yedges[ibn]],
####                    aspect='equal', interpolation='none',
####                    cmap=cmap, vmin=-0.05, vmax=0.05)
####    
####    plt.xlabel('Xpos')
####    plt.ylabel('Ypos')
####    
####    #        if(j==1):
####    plt.clim(-0.05,0.05)
####            
####    if(i==1):
####        plt.title('e1')
####    else:
####        plt.title('e2')
####        
####
####plt.gca().set_aspect('equal',adjustable='box')
####plt.tight_layout()
####        
##### Make an axis for the colorbar on the right side
####
####cax = fig.add_axes([0.1, 0.05, 0.8, 0.01])
####fig.colorbar(im, cax=cax, orientation='horizontal')
####
####plt.savefig(md+'/PLOTS/e_X_Y_%s_V1.0.0A.png' % (patch))
#####plt.show()
##
#####################################################
#####################################################
#####################################################
##
###im = plt.imshow(e1_mean, origin='lower',
###                    aspect='equal', interpolation='none',
###                    cmap=cmap) #, vmin=-0.005, vmax=0.005)
###plt.colorbar()
###plt.show()
##
###im = plt.hist(c1_map.flatten(), bins=100)
###plt.show()
##

