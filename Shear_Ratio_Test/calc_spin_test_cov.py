import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import ascii 

SPINDIR='/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/Shear_Ratio_Test/Output/'
#'/disk09/KIDS/K1000_TWO_PT_STATS//OUTSTATS/SHEAR_RATIO/'
#number of spin trials, angular bins (ntbins) and lens bins (nlbins)

JZBIN=1  # redshift bin of the sources
ntrials=400
ntbins=4
lbins=[3,4,5]

gtlens=np.zeros([ntbins,ntrials])

for ilens in range(len(lbins)):
    for ispin in range(ntrials):
        gtfile=SPINDIR+'SPIN/K1000_GT_SPIN_'+str(ispin)+'_6Z_source_'+str(JZBIN)+'_5Z_lens_'+str(lbins[ilens])+'.asc'
        gtspindat=ascii.read(gtfile)
        gtlens[:,ispin]=gtspindat['gamT'] 
    if ilens==0:      
        gtprev = np.copy(gtlens)
    else:    
        gtall = np.vstack((gtprev,gtlens))
        gtprev= np.copy(gtall)

Cov=np.cov(gtall)
Corr=np.corrcoef(gtall)

#plt.imshow(Cov, interpolation='None')
plt.imshow(Corr, interpolation='None')
plt.colorbar()
plt.axis('off')
plt.show()
plt.savefig('Shear_Ratio_Correlation_source_bin'+str(JZBIN)+'.png')
