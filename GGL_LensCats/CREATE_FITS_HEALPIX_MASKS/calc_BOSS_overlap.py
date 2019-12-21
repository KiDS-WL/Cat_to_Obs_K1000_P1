# quick-write python script to calculate the effective BOSS overlap area from the BOSS randoms
# Read in the BOSS randoms that have had the KiDS Mask appended in makelenscats.py

from astropy.io import fits
import numpy as np

import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d

cmap = plt.cm.jet
cmap.set_bad('w', 1.)

DD='/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS'

#To apply the 9-band photo-zs mask use
#bitmask=0x6FFC

#You might however just want to know where the gri information is
#for this you want
#bitmask=0x681C

#Or maybe you just want to know if the VST has looked at this area
bitmask=0x4000

#Lets first calculate how many randoms were in the original NGP random files
nran_ngp = np.zeros(2)
for iran in range(2):
    
    #you might assume that you should  use the full CMASS/LOWZ random file
    #this won't work though as this catalogue has variable density as LOWZ
    #has a different footprint from CMASSS
    #original_randomfile=DD+'/BOSS_original/random'+str(iran)+'_DR12v5_CMASSLOWZTOT_North.fits'
    
    #It's therefore best to analyse the surveys seperately
    original_randomfile=DD+'/BOSS_original/random'+str(iran)+'_DR12v5_CMASS_North.fits'
    #original_randomfile=DD+'/BOSS_original/random'+str(iran)+'_DR12v5_LOWZ_North.fits'
    
    hdulist = fits.open(original_randomfile)
    randomcat = hdulist[1].data
    zspec=randomcat.field('Z')

    #Split the random sample into the two z-bins in Sanchez et al
    nran_ngp[0] = nran_ngp[0] + len(zspec[(zspec <=0.5) & (zspec>0.2)])
    nran_ngp[1] = nran_ngp[1] + len(zspec[(zspec <= 0.75) & (zspec>0.5)])
    print (nran_ngp)

#Now lets see how many remain once the KiDS mask has been applied
for ilens in range(1):

    #randomfile=DD+'/BOSS_random_z'+str(ilens+1)+'.fits'
    randomfile=DD+'/BOSS_random_CMASS_z'+str(ilens+1)+'.fits'
    #randomfile=DD+'/BOSS_random_LOWZ_z'+str(ilens+1)+'.fits'
    
    #Read in the random catalogue and the MASK
    hdulist = fits.open(randomfile)
    randomcat = hdulist[1].data
    KIDSMASK=randomcat.field('KIDSMASK')
    RAall=randomcat.field('ALPHA_J2000')
    Decall=randomcat.field('DELTA_J2000')
    
    #filter based on the 9-band, gri mask, or wcs mask
    ifilter=np.logical_not(np.array(KIDSMASK.astype(int) & bitmask, dtype=bool))
    not_in_mask=KIDSMASK[ifilter]
    RA=RAall[ifilter]
    Dec=Decall[ifilter]

    #Generally a good idea to make a plot for sanity checks
    nbins=600
    plt.figure(figsize=(8,7))
    nran_binned, xedges, yedges, binnum = binned_statistic_2d(Dec, RA, not_in_mask,statistic='count', bins=nbins)
    plt.imshow(nran_binned, origin='lower', 
                        extent=[yedges[0], yedges[nbins],
                                xedges[0], xedges[nbins]],
                        aspect='auto', interpolation='none', cmap=cmap)
    plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('Dec')
    #plt.savefig('random_distribution_ilens_%d.png' % (ilens))
    plt.savefig('CMASS_random_distribution_ilens_%d_KiDS_wcs.png' % (ilens))
    #plt.savefig('LOWZ_random_distribution_ilens_%d.png' % (ilens))
    plt.show()

    #We now compare the number of randoms within the KIDS footprint
    #to the number in the full BOSS area catalogue in the two redshift slices

    frac_in_footprint=len(not_in_mask)*1.0 / nran_ngp[ilens]

    #What is the effective area in the NGP?
    #The total area is 9329 degrees from Alam et al (arxiv: 1607.03155): 9329 sq degrees
    #We want the NGP effective area though
    #- for this use Table 2 of Reid et al (arxiv: 1509.06529): 6851 sq degrees 

    print ('ilens %d BOSS overlap area %f'%(ilens,frac_in_footprint*6851.0))

