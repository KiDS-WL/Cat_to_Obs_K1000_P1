# 29/11/2019, B. M. Giblin
# Make n(z) source histograms from the MICE2 mocks
# in the range (0.1,0.3,0.5,0.7,0.9,1.2) and save them.
# These are to be used in the shear ratios tests.
# ALSO save separate shear catalogues for these tomo bins to be used in the SRT.

import numpy as np
import os
from astropy.io import fits
from scipy.integrate import simps

LENS_TYPE = 'GAMA_data'
numz = 5



for i in range(1, numz+1):
    f = fits.open('%s/lens_cat_%sZ_%s.fits' %(LENS_TYPE,numz,i))
    z = f[1].data['Z']
    weight = f[1].data['WEICOMP']
    f.close()
    bin_lo_edges = np.linspace( 0., 0.7, 120)
    bin_edges = np.append( bin_lo_edges, bin_lo_edges[-1]+(bin_lo_edges[1]-bin_lo_edges[0]) )
    histo_Z = np.histogram(z, bins=bin_edges, weights=weight)[0]

    # There's something weird with the GAMA nofz, huge spike at z=0: remove
    # (this might be causing the very poor p-value when we do K1000xGAMA incl. IAs).
    if LENS_TYPE == 'GAMA_data':
        histo_Z[0] = histo_Z[1]

    area = simps( histo_Z, bin_lo_edges )
    histo_Z = histo_Z / area
    # Save it as ascii or FITS?
    # FITS
    hdulist = fits.BinTableHDU.from_columns(
            [fits.Column(name='binstart', format='1D', array=bin_lo_edges),
             fits.Column(name='density', format='1D', array=histo_Z) ])
    hdulist.writeto('%s/nofz_cat_%sZ_%s.fits' %(LENS_TYPE,numz,i), overwrite=True)
    # ascii
    np.savetxt('%s/nofz_cat_%sZ_%s.dat' %(LENS_TYPE,numz,i),
               np.c_[bin_lo_edges, histo_Z], header='# z_spec_bin, n(z)' )



def Plot_Nz():
    import matplotlib.pyplot as plt
    for i in range(1, numz+1):
        tmp_bin, tmp_hist = np.loadtxt('%s/nofz_cat_%sZ_%s.dat' %(LENS_TYPE,numz,i),
                                   usecols=(0,1), unpack=True)
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
    plt.xlabel('Zspec')
    plt.ylabel('N(z)')
    plt.title('%s' %LENS_TYPE)
    plt.show()
            
    return
Plot_Nz()
