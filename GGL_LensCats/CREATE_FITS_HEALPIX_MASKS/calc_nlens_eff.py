# quick-write python script to calculate the effective number density of lenses from BOSS
# or from 2dFLenS
# Read in the BOSS randoms that have had the KiDS Mask appended in makelenscats.py

from astropy.io import fits
import numpy as np

DD='/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS'

#Choose the survey
survey='2dFLenS'
#survey='BOSS'

#Choose the number of KIDS bands for the overlap region
#nbands=9  # full 9-band area
nbands=3 # gri area
#nbands=1  # in the KiDS imaging area - i.e wcs mask

if nbands==9:
    #To apply the 9-band photo-zs mask use
    bitmask=0x6FFC

    if survey=='BOSS':
        #BOSS 9-band overlap area
        area=322.255634
    else:
        #2dFLenS 9-band overlap area
        area=342.879925
elif nbands==3:
    #You might however just want to know where the gri information is
    #for this you want
    bitmask=0x681C
    if survey=='BOSS':
        #BOSS 3-band overlap area
        area=340.132473
    else:
        #2dFLenS 3-band overlap area
        area=355.728369 
elif nbands==1:
    #You might however just want to know where the VST has looked
    #for this you want
    bitmask=0x4000
    if survey=='BOSS':
        #BOSS 1-band overlap area
        area=409.189476
    else:
        #2dFLenS 1-band overlap area
        area=425.183122

#Now lets see how many lenses remain once the KiDS mask has been applied

for ilens in range(2):
    datafile=DD+'/'+survey+'_data_z'+str(ilens+1)+'.fits'
    
    #Read in the data catalogue and the MASK
    hdulist = fits.open(datafile)
    datacat = hdulist[1].data
    KIDSMASK=datacat.field('KIDSMASK')
    Weightall=datacat.field('WEICOMP')  # BOSS and 2dFLens without gri weighting
    
    if (survey=='2dFLenS' and nbands==3):
        Weightall=datacat.field('WEIMAG') #2dFLenS with gri weighting

    #filter based on the 9-band or gri mask
    ifilter=np.logical_not(np.array(KIDSMASK.astype(int) & bitmask, dtype=bool))
    weight=Weightall[ifilter]

    #to get a neff we need to sum the weights, in the same way we do for lensfit
    #neff = (1/A)*(Sum weights)^2/(Sum weights^2)
    neff = ((np.sum(weight))**2/np.sum(weight*weight))/(area*60*60)

    print('ilens %d, %d-band, %s lenses neff per sq arcmin = %f'%(ilens, nbands, survey,neff))
