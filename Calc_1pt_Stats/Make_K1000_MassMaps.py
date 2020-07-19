# 15/07/2020, B. M. Giblin, Postdoc, Edinburgh
# Read in shear data for KiDS-1000 and perform spherical harmonic transform
# to get kappa maps. Also do a series of noise realisations to make
# signal/noise maps.

import sys
import numpy as np
import healpy as hp
from matplotlib import rcParams
from astropy.io import fits
import matplotlib.pyplot as plt
import time

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}

plt.rc('font', **font)


# Read in user input to set the patch, blind, zmin,zmax, nbootstrap

if len(sys.argv) <4:
    print( "Usage: %s Patch Blind ZBmin ZBmax" % sys.argv[0] )
    sys.exit(1)
else:
    NorS=sys.argv[1]  # N for North, S for South, or 'All' for both.
    Blind=sys.argv[2] # Just put 'A'
    ZBlo=sys.argv[3]  # If numbers are inputted, this Z_B cut is applied to data
    ZBhi=sys.argv[4]  # For no ZBcut, put 'None' for both of these 



def Pickle_Data(ra,dec,e1,e2,w,ZB):
    print( "Pickling RA,Dec,e1,e2,w,ZB...." )
    outputrd='Catalogues/MassMapFile-K%s.Blind%s.ra_dec_e1_e2_w_ZB' %(NorS,Blind)
    np.save(outputrd, np.column_stack((ra,dec,e1,e2,w,ZB)) )
    return

nside = 1024 # ps ~ 3 arcmin/pxl
ps = 60*(4*np.pi * (180./np.pi)**2 / (12*nside**2))**0.5 # pixel scale in arcmin / pxl
sigma = 10. # smoothing scale in arcmin per pixel
Read_Cat_Or_Pickle = None # "Cat", "Pickle", None    
Run_Or_Read_Noise = "Read"  # Either run the noise realisations, or read in pre-saved
Noise_start = 0
Noise_end = 100
Num_Noise = Noise_end - Noise_start

if Read_Cat_Or_Pickle == "Cat":
    inDIR = '/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/'

    f = fits.open('%s/K1000_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses.cat' %(inDIR,NorS))
    if NorS == 'ALL':
        iext = 2
    else:
        iext = 1
    ra = f[iext].data['ALPHA_J2000']
    dec = f[iext].data['DELTA_J2000']
    ZB = f[iext].data['Z_B']
    # Not certain which e1/e2/w to read in?
    # ...this?
    #mce1 = f[iext].data['raw_e1']
    #mce2 = f[iext].data['raw_e2']
    #mcm1 = 1+f[iext].data['metacal_m1']
    #mcm2 = 1+f[iext].data['metacal_m2']
    #mcw=(mcm1+mcm2)*0.5
    # ...or this?
    #autocal columns                                                                                                                   
    e1 = f[iext].data['autocal_e1_%s' %Blind]
    e2 = f[iext].data['autocal_e2_%s' %Blind]
    w = f[iext].data['recal_weight_%s' %Blind]
    Pickle_Data(ra,dec,e1,e2,w,ZB)
    
elif Read_Cat_Or_Pickle == "Pickle":
    ra,dec,e1,e2,w,ZB = np.load('Catalogues/MassMapFile-K%s.Blind%s.ra_dec_e1_e2_w_ZB.npy' %(NorS,Blind)).transpose()[[0,1,2,3,4,5],:]


# Cut redshift
if str(ZBlo) == "None":
    ZBlabel = 'ZBcutNone'
    
elif float(ZBlo)<float(ZBhi):
    ZBlabel = 'ZBcut%s-%s' %(ZBlo,ZBhi)
    if Read_Cat_Or_Pickle == "Cat" or Read_Cat_Or_Pickle == "Pickle":
        print( "Making the ZBcut in the range %s to %s" %(ZBlo, ZBhi) )
        idx = np.where( np.logical_and(ZB>float(ZBlo), ZB<float(ZBhi)) )[0]
        ra=ra[idx]
        dec=dec[idx]
        #metacal columns
        #mce1=mce1[idx]
        #mce2=mce2[idx]
        #mcw=mcw[idx]
        #mcm1=mcm1[idx]
        #mcm2=mcm2[idx]
        #autocal columns
        e1=e1[idx]
        e2=e2[idx]
        w=w[idx]
        ZB=ZB[idx]


def Make_Arrays_Into_Healpix(nside, ra_array, dec_array, weight_array, other_array):
    phi = ra_array*np.pi/180. # in rads
    theta = (90. - dec_array)*np.pi/180. # in rads
    pix = hp.pixelfunc.ang2pix(nside, theta, phi, nest=False)

    hp_array = np.zeros(12*nside**2)
    hp_weight = np.zeros(12*nside**2)
    hp_counts = np.zeros(12*nside**2) # so you can take the average weight for each healpix pixel too

    # weighted average of the original array values that share the same healpix pixel
    # also create a weight healpix map for use in masking & taking weighted AVERAGE
    for i in range(0, len(pix)):
        hp_array[pix[i]] += other_array[i]*weight_array[i]
        hp_weight[pix[i]] += weight_array[i]
        hp_counts[pix[i]] += 1

    # average
    hp_array[ np.where( hp_weight != 0.)[0] ] /= hp_weight[ np.where( hp_weight != 0.)[0] ]
    hp_weight[ np.where(hp_weight != 0.)[0] ] /= hp_counts[ np.where( hp_weight != 0.)[0] ]

    # now mask
    hp_array[ np.where( hp_weight == 0.)[0] ] = hp.UNSEEN
    hp_weight[ np.where( hp_weight == 0.)[0] ] = hp.UNSEEN
    return hp_array, hp_weight


def Make_KappaMap():
    e1_hp, weight_hp = Make_Arrays_Into_Healpix(nside, ra, dec, w, e1) 
    e2_hp, weight_hp = Make_Arrays_Into_Healpix(nside, ra, dec, w, e2)

    # Make weighted versions of the ellipticity
    e1_hp_w = e1_hp * weight_hp
    e1_hp_w[ weight_hp==hp.UNSEEN ] = hp.UNSEEN
    e2_hp_w = e2_hp	* weight_hp
    e2_hp_w[ weight_hp==hp.UNSEEN ] = hp.UNSEEN

    # Make a mass map
    print( "Making kappa map" )

    t0=time.time()
    _, Elm, Blm = hp.sphtfunc.map2alm([np.zeros_like(e1_hp_w), e1_hp_w, e2_hp_w], pol=True, use_weights=True)

    kappa_hp = hp.sphtfunc.alm2map( Elm, nside, sigma=sigma *np.pi/(180.*60.) )
    kappa_hp[ weight_hp==hp.UNSEEN ] = hp.UNSEEN
    hp.write_map('MassMaps/K%s_Kappa_nside%s_sigma%sarcmin_%s.fits' %(NorS, nside,sigma,ZBlabel),
              kappa_hp, overwrite=True)

    t1 = time.time()
    print( "Making the kappa map took %.1f s" %(t1-t0) )
    return kappa_hp
# Make Map
#kappa_hp = Make_KappaMap()
# Read pre-saved map
kappa_hp = hp.read_map('MassMaps/K%s_Kappa_nside%s_sigma%sarcmin_%s.fits' %(NorS, nside,sigma,ZBlabel))
mask = kappa_hp==hp.UNSEEN

def Plot_Mollview(Map, rot, mini, maxi, title):
    plt.figure(figsize=(9,9))
    hp.mollview(Map, fig=1, title=title , min=mini, max=maxi, cbar=False, rot=rot)
    hp.graticule(coord='G', color='white')
    plt.show()
    return

def Plot_Gnomview(Map, figsize, rot, reso, xsize, ysize, mini, maxi, title, savename):
    # gnomview seems to always have ra as straight horizontal line
    # sonot very suitable for South field
    plt.figure(figsize=figsize)
    hp.gnomview(Map, fig=1, min=mini, max=maxi, cbar=False, cmap='magma', badcolor='black',
                rot=rot, reso=ps, xsize=xsize, ysize=ysize,
                title=title, notext=True, format='%.1f' )
    #hp.graticule(coord='C', color='white')
    #plt.savefig(savename)
    plt.show()
    return

def Plot_Orthview(Map, rot, mini, maxi, savename):
    plt.figure(figsize=(12,12))
    hp.orthview(Map, fig=1, min=mini, max=maxi, cbar=True,
                rot=rot, half_sky=True,xsize=5000, cmap='inferno',
                title=None, format='%.1f', flip='geo' )
    hp.graticule(dpar=10, dmer=10, coord='G', color='white')
    plt.savefig(savename) 
    plt.show()
    return


# Gnomview
#Plot_Gnomview(kappa_hp, (20,1.5), [182.,0.,0.], ps, 3000, 200, -0.1, 0.1, 'KiDS-North',
#              'MassMaps/gnomview_KN_Kappa_nside%s_sigma%sarcmin_%s.png' %(nside,sigma,ZBlabel) )
#Plot_Gnomview(kappa_hp, (10,2.6), [12.,-36.,0.], ps, 1500, 350, -0.1, 0.1, 'KiDS-South',
#              'MassMaps/gnomview_KS_Kappa_nside%s_sigma%sarcmin_%s.png' %(nside,sigma,ZBlabel) )

# Orthview
#Plot_Orthview(kappa_hp, [0.,0.,0.], -0.1, 0.1, 'KiDS-1000')
# Mollview
#Plot_Mollview(kappa_hp, [182.,0.,0.], -0.1, 0.1, 'KiDS-1000')


# !!! EITHER CREATE OR READ IN NOISE REALISATIONS !!!
t1 = time.time()
kappa_rng_array = np.zeros([ Num_Noise, int(12*nside**2) ])
if Run_Or_Read_Noise == "Run":
    SN_level = 0.27 # avg on sigma_e across 5 bins
    k=0
    for i in range( Noise_start, Noise_end ):
        print( "Making noise realisation %s" %i )
        np.random.seed(i)
        e1_rng = np.random.normal(0., SN_level, int(12*nside**2))
        np.random.seed(i+1000)
        e2_rng = np.random.normal(0., SN_level, int(12*nside**2))
        # re-mask
        e1_rng[ mask ] = hp.UNSEEN
        e2_rng[ mask ] = hp.UNSEEN

        _, Elm_rng, _ = hp.sphtfunc.map2alm([np.zeros_like(e1_rng), e1_rng, e2_rng], pol=True, use_weights=True)
        kappa_rng = hp.sphtfunc.alm2map( Elm_rng, nside, sigma=sigma *np.pi/(180.*60.) )
        kappa_rng[ mask ] = hp.UNSEEN
        hp.write_map('MassMaps/K%s_NoiseKappa%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,i,nside,sigma,ZBlabel),
                     kappa_rng, overwrite=True)
        kappa_rng_array[k,:] = kappa_rng
        k+=1
        
    t2 = time.time()
    print( "Making %s noise maps took %.1f s" %(Num_Noise,(t2-t1)) )

    kappa_rng_std = np.std( kappa_rng_array, axis=0 )
    kappa_rng_std[ mask ] = hp.UNSEEN
    hp.write_map('MassMaps/K%s_NoiseKappaStd%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel),
                 kappa_rng_std, overwrite=True)
    
elif Run_Or_Read_Noise == "Read":
    kappa_rng_std = hp.read_map('MassMaps/K%s_NoiseKappaStd%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel))
        
else:
    print( "Not running or reading any noise realisations, so exiting code." )
    


SNR = kappa_hp / kappa_rng_std
SNR[ mask ] = hp.UNSEEN
#Plot_Gnomview(SNR, (10,3), [12.,-36.,0.], ps, 1500, 350, -1, 1, None,
#              'MassMaps/gnomview_KS_SNR%s-%s_nside%s_sigma%sarcmin_%s.png' %(Noise_start,Noise_end, nside,sigma,ZBlabel) )
Plot_Orthview(SNR, [10.,-36.,0.], -5, 5,
              'MassMaps/orthview_KS_SNR%s-%s_nside%s_sigma%sarcmin_%s.png' %(Noise_start,Noise_end, nside,sigma,ZBlabel) )
Plot_Orthview(SNR, [182.,0.,0.], -5, 5,
              'MassMaps/orthview_KN_SNR%s-%s_nside%s_sigma%sarcmin_%s.png' %(Noise_start,Noise_end, nside,sigma,ZBlabel) )
