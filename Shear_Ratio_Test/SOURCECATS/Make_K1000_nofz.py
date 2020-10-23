# 27/03/2020, B. M. Giblin
# Make n(z) source histograms for K1000
# in the range (0.1,0.3,0.5,0.7,0.9,1.2) and save them.

import numpy as np
import os
from astropy.io import fits

SOMFLAGNAME = 'Fid'
Blind = 'C'
OutlierPeak = True                # If set to true, create an outlier peak
                                   # in the tomo bins specified by OL_bins
                                   # of sizes given by OL_frac
OL_bins = '12345'                      # bins to add an outlier peak, listed in order: e.g., '123', 
OL_frac = [0.10, 0.10, 0.10, 0.10, 0.10]   # fraction of sources you want to be represented in the peak

numz = 5
inDIR = '/disk09/KIDS/K1000_TWO_PT_STATS/SOM_NofZ/'
outDIR = 'SOM_NofZ/' 


for i in range(1, numz+1):
    filename = 'K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_%s_blind%s_TOMO%s_Nz'%(SOMFLAGNAME, Blind,i)
    f = fits.open('%s/%s.fits' %(inDIR,filename))
    z = f[1].data['binstart']
    density = f[1].data['density']
    f.close()

    # ascii
    np.savetxt('%s/%s.asc' %(outDIR,filename),
               np.c_[z, density], header='# binstart, density' )

    # Also make a FITS copy of the nofz if OutlierPeak is set to True
    # this is an experiment to see if adding an outlier fraction
    # in the form of a bulge in the nofz at high-z, improves the results
    # of the SRT.
    if OutlierPeak:
        filename_ol = filename + '_OutlierPeaksInBins%s' %OL_bins
        density_wpeak = np.copy( density )
        if str(i) in OL_bins:
            print("Creating a %s outlier fraction at high-z in bin %s" %(OL_frac[i-1],i) )
            # Add a peak of the specified size at z>1.45
            idx_peak = np.where( z>1.45 )[0][0]
            sum_heights = np.sum( density[ z != z[idx_peak] ] ) # Find the sum-total density
            density_wpeak[ idx_peak ] = OL_frac[i-1] * sum_heights

        # FITS (saved for Dls_over_Ds.py to use)
        hdulist = fits.BinTableHDU.from_columns(
            [fits.Column(name='binstart', format='1D',array=z),
             fits.Column(name='density', format='1E', array=density_wpeak)])
        hdulist.writeto('%s/%s.fits' %(outDIR,filename_ol), overwrite=True)
        # ascii (saved for MakeNofZForCosmosis_LensAndSource.py in the kcap directory to read)
        np.savetxt('%s/%s.asc' %(outDIR,filename_ol),
               np.c_[z, density], header='# binstart, density' )

        
def Plot_Nz():
    import matplotlib.pyplot as plt
    for i in range(1, numz+1):
        filename = 'K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_%s_blind%s_TOMO%s_Nz'%(SOMFLAGNAME, Blind,i)
        if OutlierPeak:
            filename_ol = filename + '_OutlierPeaksInBins%s' %OL_bins
            f = fits.open( '%s/%s.fits' %(outDIR,filename_ol) )
            tmp_bin = f[1].data['binstart']
            tmp_hist = f[1].data['density']
            f.close()
            title = 'Outlier peak added at high-z in bins %s' %OL_bins
        else:
            tmp_bin, tmp_hist = np.loadtxt('%s/%s.asc' %(outDIR,filename),
                                   usecols=(0,1), unpack=True)
            title=None
            
        if i == 1:
            hist = np.reshape(tmp_hist, (1,len(tmp_hist)))
            bins = np.reshape(tmp_bin, (1,len(tmp_bin)))
        else:
            hist = np.vstack((hist, tmp_hist))
            bins = np.vstack((bins, tmp_bin))
        #print( i, hist.shape )

    plt.figure()
    colors = ['darkblue', 'orange', 'magenta', 'lawngreen', 'dimgrey']
    for i in range(hist.shape[0]):
        plt.bar( bins[i], hist[i], width=(bins[i,1]-bins[i,0]), color=colors[i], edgecolor=colors[i], alpha=0.5 )
    plt.xlabel('Z')
    plt.ylabel('N(z)')
    plt.title(title)
    plt.show()
    return
Plot_Nz()
