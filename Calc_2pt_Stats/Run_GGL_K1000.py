# 30/08/2019: B. M. Giblin, PhD Student, Edinburgh
# Code to calculate the GGL signal about BOSS gals with K1000
# Do this in quadrants to test if there is residual additive shear bias.

import numpy as np
import treecorr
import sys
from astropy.io import fits
import os
import time

# Plotting modules
#from matplotlib import rc
import matplotlib as plt
#rc('text',usetex=True)
#rc('font',size=24)
#rc('legend',**{'fontsize':30})
#rc('font',**{'family':'serif','serif':['Computer Modern']})

bin_type = "Linear"		# If "Linear" calculates the normal gamma_T(theta) where theta is log-spaced
                                # If "TwoD" calculates gamma_T (delta_x,delta_y) about lenses at (0,0).

# INPUTS FROM THE COMMAND LINE
NorS = sys.argv[1]    # N for North, S for South
Blind = sys.argv[2]   # Just put 'A' for now

ZBlo_s = sys.argv[3]  # If numbers are inputted, this Z_B cut is applied to sources
ZBhi_s = sys.argv[4]  # For no ZBcut, put None for both of these

Zlo_l = sys.argv[5]  # If numbers are inputted, this Z cut is applied to lenses
Zhi_l = sys.argv[6]  # For no Zcut, put None for both of these

print("Reading in the source data...")
inDIR = '/home/bengib/KiDS1000_NullTests/'
data = np.load('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy'%(inDIR,NorS,Blind))
ra_s,dec_s,e1_s,e2_s,w_s,ZB_s = data.transpose()[[0,1,2,3,4,5],:]

print("Making a redshift cut on the source data, if specified to do so.")
if str(ZBlo_s) == "None":
    ZBlabel = "None"     # NB: ZBlabel applies to sources, Zlabel to lenses
else:
    ZBlabel = "%s-%s" %(ZBlo_s, ZBhi_s)
    idx = np.where( np.logical_and(ZB_s>float(ZBlo_s), ZB_s<float(ZBhi_s)) )[0]
    ra_s=ra_s[idx]
    dec_s=dec_s[idx]
    e1_s=e1_s[idx]
    e2_s=e2_s[idx]
    w_s=w_s[idx]
    ZB_s=ZB_s[idx]

def Read_In_Lens_Data(fitsname):
    f = fits.open(fitsname)
    ra = f[1].data['ALPHA_J2000']
    dec = f[1].data['DELTA_J2000']
    Z = f[1].data['Z']
    w_comp = f[1].data['WEICOMP ']
    w_fkp = f[1].data['WEICOMP ']
    w = w_comp*w_fkp
    return ra, dec, w, Z

z_l = '1'	# Redshift of lenses
print("Reading in the lens data...")
ra_l, dec_l, w_l, Z_l = Read_In_Lens_Data('/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/BOSS_data_z%s.fits' %z_l)
print("Reading in the randoms data...")
ra_r, dec_r, w_r, Z_r = Read_In_Lens_Data('/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/BOSS_random_z%s.fits' %z_l)

print("Making a redshift cut on the lens data, if specified to do so.")
if str(Zlo_l) == "None":
    Zlabel = "None"   # NB: ZBlabel applies to sources, Zlabel to lenses
else:
    Zlabel = "%s-%s" %(Zlo_l, Zhi_l)
    # Lenses
    idx = np.where( np.logical_and(Z_l>float(Zlo_l), Z_l<float(Zhi_l)) )[0]
    ra_l=ra_l[idx]
    dec_l=dec_l[idx]
    w_l = w_l[idx]
    Z_l = Z_l[idx]
    # Randoms
    idx = np.where( np.logical_and(Z_r>float(Zlo_l), Z_r<float(Zhi_l)) )[0]
    ra_r=ra_r[idx]
    dec_r=dec_r[idx]
    w_r = w_r[idx]
    Z_r = Z_r[idx]


def Run_GGL(cat1, cat2):
    if bin_type == "Linear":
        # Need to provide the sep_units as arcmin
        ng = treecorr.NGCorrelation(bin_type=bin_type, min_sep=min_sep, max_sep=max_sep,
                                    bin_slop=bin_slop, nbins=nbins, sep_units=sep_units,
                                    metric=metric)
    elif bin_type == "TwoD":
        ng = treecorr.NGCorrelation(bin_type=bin_type, min_sep=min_sep, max_sep=max_sep,
                                    bin_slop=bin_slop, nbins=nbins, 
                                    metric=metric)
    print("Calculating GGL signal...")
    ng.process(cat1, cat2)
    return ng.rnom, ng.meanr, ng.xi, ng.xi_im, ng.weight


def Filter_Patches(idx, Patches_tmp, ra_lo_tmp, ra_hi_tmp, dec_lo_tmp, dec_hi_tmp):
    Patches_tmp = np.array(Patches_tmp)[idx]
    ra_lo_tmp = ra_lo_tmp[idx]
    ra_hi_tmp = ra_hi_tmp[idx]
    dec_lo_tmp = dec_lo_tmp[idx]
    dec_hi_tmp = dec_hi_tmp[idx]
    return  Patches_tmp, ra_lo_tmp, ra_hi_tmp, dec_lo_tmp, dec_hi_tmp


# gamma_t falls off at about 20 arcsec from lens
# Measure out to 20 arcmin to be sure.
# Will see freaky PSF effects down to 0 pxls --> that's okay.
min_sep_arcmin = 0.
max_sep_arcmin = 20


print("Setting up the catalogues")
if bin_type == "Linear":
    # GIBLIN YOU WANT THIS TO MATCH 2D BINNING AS MUCH AS POSS.
    # THEREFORE MATCH ANG.RANGE, DON'T DO LOG-SPACED.
    # 1D log-spaced angular sep. bins
    nbins=21
    min_sep=min_sep_arcmin   # arcmin
    max_sep=2*max_sep_arcmin  # arcmin
    sep_units = 'arcmin' 

    metric='Arc'
    bin_slop=0.05 # 0.1/(np.log(max_sep/min_sep)/float(nbins))

    # Source catalogue
    cat_s = treecorr.Catalog(ra=ra_s, dec=dec_s,
                        ra_units='degrees', dec_units='degrees',
                        g1=e1_s, g2=e2_s, w=w_s)
    # Lens catalogue
    cat_l = treecorr.Catalog(ra=ra_l, dec=dec_l,
                        ra_units='degrees', dec_units='degrees',
                        w=w_l)
    # Randoms catalogue
    cat_r = treecorr.Catalog(ra=ra_r, dec=dec_r,
                        ra_units='degrees', dec_units='degrees',
                        w=w_r)
    theta, mean_theta, xi, xi_im, weight_xi = Run_GGL(cat_l, cat_s)
    theta_r, mean_theta_r, xi_r, xi_im_r, weight_xi_r = Run_GGL(cat_l, cat_s)    
    np.savetxt('Results/K%s.Blind%s.gamma_tx.SourceZBcut%s.LensZcut%s.dat'%(NorS,Blind,ZBlabel,Zlabel),
               np.c_[theta, mean_theta, xi, xi_im, weight_xi],
               header='# theta [arcmin], mean theta [arcmin], gamma_t, gamma_x, weight')
    np.savetxt('Results/K%s.Blind%s.gamma_tx_randoms.SourceZBcut%s.LensZcut%s.dat'%(NorS,Blind,ZBlabel,Zlabel),
               np.c_[theta_r, mean_theta_r, xi_r, xi_im_r, weight_xi_r],
               header='# theta [arcmin], mean theta [arcmin], gamma_t[about randoms], gamma_x[about randoms], weight')
        
        
elif bin_type == "TwoD":
	# 2D lin-spaced delta_x,delta_y bins
	# See https://rmjarvis.github.io/TreeCorr/_build/html/binning.html#twod
    from Functions_2DGGL import Select_Patch
    from astropy.io import fits
    from astropy.wcs import WCS
    nbins = 21
    metric = 'Euclidean'
    # pxl scale in fits file is 0.213 arcsec / pxl,
    min_sep = min_sep_arcmin * 60 / 0.213
    max_sep = max_sep_arcmin * 60 / 0.213   # The maximum absolute value of dx or dy to include
    # bin_size = # The width of the bins in dx and dy
    bin_slop=0.05 # Recommended by error message

    # Load the patch information to split BOSS galaxies up.
    Patch_File = 'KiDS_Patches.dat'
    Patches = []
    f = open(Patch_File)
    lines = f.readlines()
    for x in lines:
        Patches.append(x.split(' ')[0])
    f.close()
    ra_lo, ra_hi, dec_lo, dec_hi = np.loadtxt(Patch_File, usecols=(1,2,3,4), unpack=True)

    # Forget about all the South patches (cuts 1510 -> 782)
    idx_N = np.where(dec_lo > -10)[0]
    Patches, ra_lo, ra_hi, dec_lo, dec_hi = Filter_Patches(idx_N, Patches, ra_lo, ra_hi, dec_lo, dec_hi)
    # Also some of these Patches don't have corresponding observed tiles
    # (cuts 782 -> 519)
    idx_exists = []
    for i in range(len(Patches)):
        if os.path.exists('/disk09/KIDS/KIDSCOLLAB_V1.0.0/%s'%Patches[i]):
            idx_exists.append(i)
    Patches, ra_lo, ra_hi, dec_lo, dec_hi = Filter_Patches(idx_exists, Patches, ra_lo, ra_hi, dec_lo, dec_hi)
        

    Ncats = len(Patches)
    Nlenses = np.zeros(Ncats)                      # Keep track of lenses per patch
    Nsources = np.zeros_like(Nlenses)              # Keep track of sources per patch
    Nrandoms = np.zeros_like(Nlenses)              # Keep track of randoms per patch
    idx_goodtiles = []                             # Store which tiles contain lens+sources
    idx_r_goodtiles = []                           # Store which tiles contain randoms+sources
    
    # Arrays for xi about actual lenses
    xi_store = np.empty([ Ncats, nbins, nbins ])   # Store the 2D gamma_t's
    xi_im_store = np.empty_like(xi_store)          # Stor the 2D gamma_x's
    weight_store = np.empty_like(xi_store)         # Same for the weights
    # Arrays for xi about randoms
    xi_r_store = np.empty_like(xi_store)                # Same as xi_store but with randoms
    xi_rim_store = np.empty_like(xi_store)              # Same as xi_im_store but with randoms
    weight_r_store = np.empty_like(xi_store)            # Same as weight_store but with randoms
    
    #dxdy_store = np.empty_like(xi_store)          # Checked dxdy same for every tile, no need to store.
                                                   # This contains 2D distances to lens
    dxdy_mean_store = np.empty_like(xi_store)      # store dxdy mean-centre, not used for anything.
    
    t1 = time.time()
    for i in range(Ncats):
        print("------------------------------------------------------------------------")
        print( "On Patch %s of %s"%(i,Ncats) )

        # Get the lenses in this patch.
        w_lp, ra_lp, dec_lp = Select_Patch(w_l, ra_l, dec_l,
                                       ra_lo[i],ra_hi[i],dec_lo[i],dec_hi[i])
        Nlenses[i] = len(w_lp)
        
        # Get the randoms in this patch.
        w_rp, ra_rp, dec_rp = Select_Patch(w_r, ra_r, dec_r,
                                       ra_lo[i],ra_hi[i],dec_lo[i],dec_hi[i])
        Nrandoms[i] = len(w_rp)
        
        # Get the sources in the patch.
        e1_sp, ra_sp, dec_sp = Select_Patch(e1_s, ra_s, dec_s,
                                       ra_lo[i],ra_hi[i],dec_lo[i],dec_hi[i])
        e2_sp, ra_sp, dec_sp = Select_Patch(e2_s, ra_s, dec_s,
                                       ra_lo[i],ra_hi[i],dec_lo[i],dec_hi[i])
        w_sp, ra_sp, dec_sp = Select_Patch(w_s, ra_s, dec_s,
                                       ra_lo[i],ra_hi[i],dec_lo[i],dec_hi[i])
        Nsources[i] = len(ra_sp)
        
        # Convert ra_lp,dec_lp to X,Y
        Tile_Name = '/disk09/KIDS/KIDSCOLLAB_V1.0.0/%s/r_SDSS/coadd_V1.0.0A/%s_r_SDSS.V1.0.0A.swarp.cut.fits' %(Patches[i], Patches[i])
        f = fits.open(Tile_Name)
        wcs = WCS(f[0].header)
        # Checked the following works accurately.
        X_lp,Y_lp = wcs.wcs_world2pix(ra_lp,dec_lp, 1) # The X,Y pos of lens on patch.
        X_rp,Y_rp = wcs.wcs_world2pix(ra_rp,dec_rp, 1) # The X,Y pos of randoms on patch.
        X_sp,Y_sp = wcs.wcs_world2pix(ra_sp,dec_sp, 1) # The X,Y pos of source on patch.
        f.close()


        if len(X_lp)>0 and len(X_sp)>0:  # Only attempt GGL if there's lenses/sources
            idx_goodtiles.append(i)
            # Set up the catalogues
            # Source catalogue
            cat_s = treecorr.Catalog(x=X_sp, y=Y_sp,
                                 g1=e1_sp, g2=e2_sp, w=w_sp)
            # Lens catalogue
            cat_l = treecorr.Catalog(x=X_lp, y=Y_lp,
                                 w=w_lp)
            # Randoms catalogue
            cat_r = treecorr.Catalog(x=X_rp, y=Y_rp,
                                 w=w_rp)
            dxdy, dxdy_mean_store[i,:,:], xi_store[i,:,:], xi_im_store[i,:,:], weight_store[i,:,:] = Run_GGL(cat_l, cat_s)
            
        else:
            print("Either no lenses or sources for Patch %s, number %s"%(Patches[i],i))

            
        if len(X_rp)>0 and len(X_sp)>0:  # Only attempt GGL if there's randoms/sources
            idx_r_goodtiles.append(i)
            # Set up the catalogues
            # Source catalogue
            cat_s = treecorr.Catalog(x=X_sp, y=Y_sp,
                                 g1=e1_sp, g2=e2_sp, w=w_sp)
            # Randoms catalogue
            cat_r = treecorr.Catalog(x=X_rp, y=Y_rp,
                                 w=w_rp)
            dummy_dxdy, dummy_meandxdy, xi_r_store[i,:,:], xi_rim_store[i,:,:], weight_r_store[i,:,:] = Run_GGL(cat_r, cat_s)
        else:
            print("Either no randoms or sources for Patch %s, number %s"%(Patches[i],i))
            
        print("------------------------------------------------------------------------")
        
        
    t2 = time.time()
    print("Time taken to calculate gamma_t for %s patches is %.1f s" %(Ncats,t2-t1))

    # Weighted average of the tiles
    # For actual lenses
    xi_mean = np.average( xi_store[idx_goodtiles,:,:], axis=0,
                          weights=weight_store[idx_goodtiles,:,:])
    xi_im_mean = np.average( xi_im_store[idx_goodtiles,:,:], axis=0,
                          weights=weight_store[idx_goodtiles,:,:])
    # For randoms
    xi_r_mean = np.average( xi_r_store[idx_r_goodtiles,:,:], axis=0,
                          weights=weight_r_store[idx_r_goodtiles,:,:])
    xi_rim_mean = np.average( xi_rim_store[idx_r_goodtiles,:,:], axis=0,
                          weights=weight_r_store[idx_r_goodtiles,:,:])

    
    # Save the separation, 2D gamma_t & gamma_x tiles
    # For lenses
    np.save('Results/K%s.Blind%s.2Dgamma_t.Ntiles%s.SourceZBcut%s.LensZcut%s'%(NorS,Blind,len(idx_goodtiles),ZBlabel,Zlabel),
            xi_mean)
    np.save('Results/K%s.Blind%s.2Dgamma_x.Ntiles%s.SourceZBcut%s.LensZcut%s'%(NorS,Blind,len(idx_goodtiles),ZBlabel,Zlabel),
            xi_im_mean)

    # For randoms
    np.save('Results/K%s.Blind%s.2Dgamma_t_randoms.Ntiles%s.SourceZBcut%s.LensZcut%s'%(NorS,Blind,len(idx_goodtiles),ZBlabel,Zlabel),
            xi_r_mean)
    np.save('Results/K%s.Blind%s.2Dgamma_x_randoms.Ntiles%s.SourceZBcut%s.LensZcut%s'%(NorS,Blind,len(idx_goodtiles),ZBlabel,Zlabel),
            xi_rim_mean)

    # Separation grid
    np.save('Results/K%s.Blind%s.2Ddxdy.Ntiles%s.SourceZBcut%s.LensZcut%s'%(NorS,Blind,len(idx_goodtiles),ZBlabel,Zlabel),
            dxdy)
    
    #cmap = plt.cm.rainbow
    #cmap.set_bad('w', 1.)
    #plt.figure()
    #plt.imshow(xi_mean, cmap=cmap, interpolation='none', vmin=xi_mean.min(), vmax=xi_mean.max(), origin='lower')
    #plt.colorbar(label=r'$\gamma_{\rm t}$')
    #plt.savefig('Results/xi_mean.png')
    #plt.show()







