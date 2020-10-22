# quick-write python script to calculate, plot and write out the n(z) for the BOSS and 2dFLenS samples
# we can compare the BOSS n(z) for the full NGP with the samples within the KiDS footprint
# CH:  12th Dec 2019

from astropy.io import fits
import numpy as np
from matplotlib import rcParams

import matplotlib.pyplot as plt

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
    'weight' : 'normal',
        'size'   : 19}

plt.rc('font', **font)
plt.figure(figsize=(8,7))

#This is where the catalogues live on cuillin
DD='/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS'

#To apply the 9-band photo-zs mask use
#bitmask=0x6FFC

#You might however just want to know where the gri information is
#for this you want
bitmask=0x681C

#what resolution do you want the n(z) binned with?
dz=0.01

#Now lets see the n(z) when the KiDS mask has been applied
for ilens in range(2):
    
    #input data
    bossfile=DD+'/BOSS_data_z'+str(ilens+1)+'.fits'
    twodffile=DD+'/2dFLenS_data_z'+str(ilens+1)+'.fits'
    
    #output ascii n(z)
    bossnzfile=DD+'/N_of_Z/BOSS_n_of_z'+str(ilens+1)+'_res_'+str(dz)+'.txt'
    twodfnzfile=DD+'/N_of_Z/2dFLenS_n_of_z'+str(ilens+1)+'_res_'+str(dz)+'.txt'
    twodf_w_nzfile=DD+'/N_of_Z/2dFLenS_weighted_n_of_z'+str(ilens+1)+'_res_'+str(dz)+'.txt'
    allnzfile=DD+'/N_of_Z/BOSS_and_2dFLenS_n_of_z'+str(ilens+1)+'_res_'+str(dz)+'.txt'
    
    #set up the z-range to bin over and the number of bins
    if ilens==0:
        zmin=0.2
        zmax=0.5
    else:
        zmin=0.5
        zmax=0.75
    nbins=np.int((zmax-zmin)/dz)

    #Read in the BOSS catalogue weights and the MASK
    hdulist = fits.open(bossfile)
    bosscat = hdulist[1].data
    BKIDSMASK=bosscat.field('KIDSMASK')
    bossz=bosscat.field('Z')
    bossweight = bosscat.field('WEICOMP')
    
    #filter based on the 9-band or gri mask
    ibfilter=np.logical_not(np.array(BKIDSMASK.astype(int) & bitmask, dtype=bool))

    #histogram the redshifts within the KiDS footprint
    if ilens==0:
        mylabel='BOSS in KiDS'
    else:
        mylabel=None
    n, bins, patches = plt.hist(bossz[ibfilter], nbins, normed=True, weights=bossweight[ibfilter], color='red',histtype=u'step',label=mylabel,linewidth=3)

    #write out the mean redshift
    print ('BOSS %d %s'%(ilens,np.average(bossz[ibfilter],weights=bossweight[ibfilter])))
    
    #and write out to file (note reporting the left corner of the bin here)
    np.savetxt(bossnzfile,np.c_[bins[0:nbins],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')

    #Read in the 2dFLenS catalogue and the MASK
    hdulist = fits.open(twodffile)
    twodfcat = hdulist[1].data
    TKIDSMASK=twodfcat.field('KIDSMASK')
    twodfz=twodfcat.field('Z')
    twodfweight = twodfcat.field('WEIMAG')
    twodfweightcomp = twodfcat.field('WEICOMP')
    
    #filter based on the 9-band or gri mask
    itfilter=np.logical_not(np.array(TKIDSMASK.astype(int) & bitmask, dtype=bool))

    #this with no weights

    if ilens==0:
        mylabel='2dFLenS in KiDS'
    else:
        mylabel=None
    n, bins, patches = plt.hist(twodfz[itfilter], nbins, normed=True, color='green', histtype=u'step',label=mylabel,linewidth=2.5)
    np.savetxt(twodfnzfile,np.c_[bins[0:nbins],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')

    #write out the mean redshift
    print ('2dFLenS %d %s'%(ilens,np.average(twodfz[itfilter])))
    
    #here we apply the gri weights

    if ilens==0:
        mylabel='2dFLenS weighted'
    else:
        mylabel=None

    n, bins, patches = plt.hist(twodfz[itfilter], nbins, normed=True, weights=twodfweight[itfilter], color='blue', histtype=u'step',label=mylabel,linewidth=2.5)
    np.savetxt(twodf_w_nzfile,np.c_[bins[0:nbins],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')

    #what does a combined unweighted 2dFLenS and BOSS in KiDS n(z) look like?
    allinkids = np.append(twodfz[itfilter],bossz[ibfilter])
    allweight = np.append(twodfweightcomp[itfilter],bossweight[ibfilter])


    if ilens==0:
        mylabel='BOSS and 2dFLenS'
    else:
        mylabel=None

    n, bins, patches = plt.hist(allinkids, nbins, normed=True, weights=allweight, color='orange', histtype=u'step',label=mylabel,linewidth=2.5)
    np.savetxt(allnzfile,np.c_[bins[0:nbins],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')

    #write out the mean redshift
    print ('All %d %s'%(ilens,np.average(allinkids,weights=allweight)))
    
#Lets overplot the n(z) in the original NGP data files

original_datafile=DD+'/BOSS_original/galaxy_DR12v5_CMASSLOWZTOT_North.fits'
hdulist = fits.open(original_datafile)
datacat = hdulist[1].data
zspec=datacat.field('Z')

weicp = datacat.field('WEIGHT_CP')
weinoz = datacat.field('WEIGHT_NOZ')
weisys = datacat.field('WEIGHT_SYSTOT')
weicompboss = weisys*(weinoz+weicp-1.)

zbin1filt=((zspec <=0.5) & (zspec>0.2))
nbins1=np.int((0.5-0.2)/dz)
bossallnzfile1=DD+'/N_of_Z/BOSS_NGP_n_of_z1.txt'
zbin2filt=((zspec <=0.75) & (zspec>0.5))
nbins2=np.int((0.75-0.5)/dz)
bossallnzfile2=DD+'/N_of_Z/BOSS_NGP_n_of_z2.txt'


n, bins, patches = plt.hist(zspec[zbin1filt], nbins1, normed=True, weights=weicompboss[zbin1filt], color='black', alpha=0.75,label='BOSS all',histtype=u'step',linewidth=4)
np.savetxt(bossallnzfile1,np.c_[bins[0:nbins1],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')

n, bins, patches = plt.hist(zspec[zbin2filt], nbins2, normed=True, weights=weicompboss[zbin2filt], color='black', alpha=0.75, histtype=u'step',linewidth=4)
np.savetxt(bossallnzfile2,np.c_[bins[0:nbins2],n],fmt=['%.3f','%.3f'],header='z_bin_left n_of_z')




plt.xlim(0.15,0.8)
plt.xlabel('z')
plt.ylabel('n(z)')
plt.legend(loc = 'upper left', fontsize=14)
plt.savefig('BOSS-2dFLenS-nofz.png')

plt.show()
