import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import ascii
import sys

# ----- Load Input ----- #
from Get_Input import Get_Input
paramfile = sys.argv[1]   # e.g. params_KiDSdata.dat
GI = Get_Input(paramfile)
SOURCE_TYPE = GI.Source_Type()
LENS_TYPE = GI.Lens_Type()
RANDOM_TYPE = GI.Random_Type()
SPINDIR='Output/SOURCE_'+str(SOURCE_TYPE)+'_LENS_'+str(LENS_TYPE)
# '/disk09/KIDS/K1000_TWO_PT_STATS//OUTSTATS/SHEAR_RATIO/'
# number of spin trials, angular bins (ntbins) and lens bins (nlbins)

JZBIN=1  # redshift bin of the sources
ntrials=70
ntbins=4
lbins=[1,2,3]

gtlens=np.zeros([ntbins,ntrials])

for ilens in range(len(lbins)):
    for ispin in range(ntrials):
        gtfile=SPINDIR+'/SPIN/GT_SPIN_'+str(ispin)+'_6Z_source_'+str(JZBIN)+'_3Z_lens_'+str(lbins[ilens])+'.asc'

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
plt.savefig(SPINDIR+'/SPIN/Shear_Ratio_Correlation_source_bin'+str(JZBIN)+'.png')

