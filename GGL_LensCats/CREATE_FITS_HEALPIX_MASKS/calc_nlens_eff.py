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
nbands=9
#nbands=4

if nbands==9:
    #To apply the 9-band photo-zs mask use
    bitmask=0x6FFC

    if survey=='BOSS':
        #BOSS 9-band overlap area
        area=319.506
    else:
        #2dFLenS 9-band overlap area
        area=341.888
elif nbands==4:
    #You might however just want to know where the gri information is
    #for this you want
    bitmask=0x681C
    if survey=='BOSS':
        #BOSS 4-band overlap area
        area=339.298
    else:
        #2dFLenS 4-band ovelrap area
        area=355.2831

#Now lets see how many lenses remain once the KiDS mask has been applied

for ilens in range(2):
    datafile=DD+'/'+survey+'_data_z'+str(ilens+1)+'.fits'
    
    #Read in the data catalogue and the MASK
    hdulist = fits.open(datafile)
    datacat = hdulist[1].data
    KIDSMASK=datacat.field('KIDSMASK')
    Weightall=datacat.field('WEICOMP')  # BOSS and 2dFLens without gri weighting
    
    if (survey=='2dFLenS' and nbands==4):
        Weightall=datacat.field('WEIMAG') #2dFLenS with gri weighting

    #filter based on the 9-band or gri mask
    ifilter=np.logical_not(np.array(KIDSMASK.astype(int) & bitmask, dtype=bool))
    weight=Weightall[ifilter]

    #to get a neff we need to sum the weights, in the same way we do for lensfit
    #neff = (1/A)*(Sum weights)^2/(Sum weights^2)
    neff = ((np.sum(weight))**2/np.sum(weight*weight))/(area*60*60)

    print('ilens %d, %d-band, %s lenses neff per sq arcmin = %f'%(ilens, nbands, survey,neff))
