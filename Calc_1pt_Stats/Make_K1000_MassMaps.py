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
#rcParams['ps.useafm'] = True
#rcParams['pdf.use14corefonts'] = True

#font = {'family' : 'serif',
#        'weight' : 'normal',
#        'size'   : 19}

#plt.rc('font', **font)


# Read in user input to set the patch, blind, zmin,zmax, nbootstrap

if len(sys.argv) <4:
    print( "Usage: %s Patch Blind ZBmin ZBmax" % sys.argv[0] )
    sys.exit(1)
else:
    NorS=sys.argv[1]  # N for North, S for South, or 'All' for both.
    Blind=sys.argv[2] # Just put 'A'
    ZBlo=sys.argv[3]  # If numbers are inputted, this Z_B cut is applied to data
    ZBhi=sys.argv[4]  # For no ZBcut, put 'None' for both of these 
    Noise_start=int(sys.argv[5])
    Noise_end=int(sys.argv[6])


def Pickle_Data(ra,dec,e1,e2,w,ZB):
    print( "Pickling RA,Dec,e1,e2,w,ZB...." )
    outputrd='Catalogues/MassMapFile-K%s.Blind%s.ra_dec_e1_e2_w_ZB' %(NorS,Blind)
    np.save(outputrd, np.column_stack((ra,dec,e1,e2,w,ZB)) )
    return

nside = 1024  #1024 # ps ~ 3 arcmin/pxl
ps = 60*(4*np.pi * (180./np.pi)**2 / (12*nside**2))**0.5 # pixel scale in arcmin / pxl
sigma = 10. #10. # smoothing scale in arcmin per pixel

Read_Cat_Or_Pickle = None            # "Cat", "Pickle", None
                                     # If your kappa map is already made, set to None
                                     # Else set to "Cat" to read from fits catalogues
                                     # or "Pickle" if data is already pickled (this is faster).
                                     # Advise to do Cat 1st time, then Pickle/None
                                     
Run_Or_Read_Noise = None             # "Read", "Read_Multiple", "Run", None
                                     # Either run the noise realisations,
                                     # read in the multiple noise maps already saved,
                                     # read in the single STDEV across the noise maps already saved
Read_SNR = False                      # or don't do anything and just read in the pre-saved SNR map

#Noise_start = 0
#Noise_end = 100
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
    print( "Making the kappa map took %.0f s with nside %s" %((t1-t0),nside) )
    return kappa_hp


if Read_Cat_Or_Pickle == None:
    # Read pre-saved map
    kappa_hp = hp.read_map('MassMaps/K%s_Kappa_nside%s_sigma%sarcmin_%s.fits' %(NorS, nside,sigma,ZBlabel))
else:
    # Make Map
    kappa_hp = Make_KappaMap()

mask = kappa_hp==hp.UNSEEN


def Fix_Holes(Map):
    t1 = time.time()
    print("Filling in tiny holes in survey footprint which are due to pixelisation. This takes ~70s with nside4096.")
    all_UNSEEN = np.where( Map == hp.UNSEEN )[0]
    # Cut this array down to just the N&S footprints
    # (between a 50% and 95% reduction in size and thus computation time for nside 1024, 4096 respectively).
    ra_UNSEEN, dec_UNSEEN = hp.pixelfunc.pix2ang(nside, all_UNSEEN, nest=False, lonlat=True)
    # North patch
    idx_N1 = np.where( np.logical_and(dec_UNSEEN>-10, dec_UNSEEN<10) )[0]
    idx_N2 = np.where( np.logical_and(ra_UNSEEN>120, ra_UNSEEN<240) )[0]
    idx_N = np.intersect1d( idx_N1, idx_N2 ) # elements of all_UNSEEN in the N patch
    # this next line can just be used to sanity check the pxls you've found are in the N patch:
    #ra_UNSEEN_N, dec_UNSEEN_N = hp.pixelfunc.pix2ang(nside, all_UNSEEN[idx_N], nest=False, lonlat=True)

    # South patch
    idx_S1 = np.where( np.logical_and(dec_UNSEEN>-38, dec_UNSEEN<-25) )[0]
    idx_S2 = np.where( np.logical_or(ra_UNSEEN>327, ra_UNSEEN<57) )[0]
    idx_S = np.intersect1d( idx_S1, idx_S2 )

    # cut all masked pxls to just those in the N&S patches 
    all_UNSEEN = all_UNSEEN[ np.append(idx_N,idx_S) ]
    # get all 8 neighbours for each pxl
    all_UNSEEN_neighbours = hp.pixelfunc.get_all_neighbours(nside, all_UNSEEN, nest=False)

    new_Map = np.copy(Map) # need new map to alter the pxls of
                           # whilst using old map to assess number of masked pairs.
    if nside == 4096:
        max_allowed_masked_pxls = 8
    else:
        max_allowed_masked_pxls = 4
        # number above is the boundary at which we say "there's too many masked adjacent pxls to
        # justify colouring this one in". So higher number means we do more colouring in.                                               
        # MAX ALLOWED VALUE IS 8. Any more than that and it will colour the entire map in.
        # this value was basically set empirically to improve map asthetic. 
        
    for i in range(len(all_UNSEEN)):
        # find out how many neighbouring pixels are SEEN & UNSEEN
        tmp_UNSEEN = np.where( Map[ all_UNSEEN_neighbours[:,i] ] == hp.UNSEEN )[0]
        tmp_SEEN = np.where( Map[ all_UNSEEN_neighbours[:,i] ] != hp.UNSEEN )[0]

        if len(tmp_UNSEEN)<max_allowed_masked_pxls: 
            new_Map[ all_UNSEEN[i] ] = Map[ all_UNSEEN_neighbours[tmp_SEEN[0],i] ]

    t2 = time.time()
    print( "Took %.0f sec to fill the holes with nside %s" %((t2-t1),nside) )
    return new_Map



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









def Plot_Orthview(figsize, Map, rot, mini, maxi, savename, xlimits, ylimits, ras, decs_lo, decs_up):
    fontsize=14
    plt.figure(figsize=figsize)
    hp.orthview(Map, fig=1, min=mini, max=maxi, cbar=None,
                rot=rot, xsize=2*(ras[-1] - ras[0])*60 / ps, # xsize is v important for determining pixelisation
                half_sky=True,cmap='inferno',badcolor='black',
                title=None) #, format='%.1f', flip='geo' )
    # Set x & y limits
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    
    # ras & decs_up/decs_low are the tick labels to be plotted with:
    #for r in ras: # plot the ra labels
        #hp.visufunc.projtext(r,decs_up-5,'%.0f'%r+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    # (might have to play around with ra positioning of dec labels)
    #hp.visufunc.projtext(ras[0]-10,decs_lo,'%.0f'%decs_lo+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    #hp.visufunc.projtext(ras[0]-10,decs_lo+10,'%.0f'%(decs_lo+10)+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    #hp.visufunc.projtext(ras[0]-10,decs_up,'%.0f'%decs_up+r'$^{\circ}$', lonlat=True, coord='G', color='white', fontsize=fontsize)
    #hp.visufunc.projplot(ras,decs_up,'wo',lonlat=True)

    # Plots ra/dec contours every dpar intervals in ra, dmer intervals in dec
    hp.graticule(dpar=10, dmer=10, coord='G', color='white')
    plt.savefig(savename) 
    plt.show()
    return


# Gnomview
#Plot_Gnomview(kappa_hp, (20,1.5), [182.,0.,0.], ps, 3000, 200, -0.1, 0.1, 'KiDS-North',
#              'MassMaps/gnomview_KN_Kappa_nside%s_sigma%sarcmin_%s.png' %(nside,sigma,ZBlabel) )
#Plot_Gnomview(kappa_hp, (10,2.6), [12.,-36.,0.], ps, 1500, 350, -0.1, 0.1, 'KiDS-South',
#              'MassMaps/gnomview_KS_Kappa_nside%s_sigma%sarcmin_%s.png' %(nside,sigma,ZBlabel) )

# Mollview
#Plot_Mollview(kappa_hp, [182.,0.,0.], -0.1, 0.1, 'KiDS-1000')

#ras_S = np.arange(-20,70,10)
#Plot_Orthview((12,4), kappa_hp, [10.,-30.], -.1, .1,
#              'MassMaps/orthview_KS_Kappa_nside%s_sigma%sarcmin_%s.png' %(nside,sigma,ZBlabel),
#              [-0.65,0.65], [-0.3,0.25], ras_S, -40, -20)




# !!! EITHER CREATE OR READ IN NOISE REALISATIONS !!!
t1 = time.time()

if Run_Or_Read_Noise == "Run":
    kappa_rng_array = np.zeros([ Num_Noise, int(12*nside**2) ])
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
    #del kappa_rng_array
    kappa_rng_std[ mask ] = hp.UNSEEN
    hp.write_map('MassMaps/K%s_NoiseKappaStd%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel),
                 kappa_rng_std, overwrite=True)


elif Run_Or_Read_Noise == "Read_Multiple":
    # USE THIS IF THE NOISE MAPS ARE SAVED BUT YOU NEED TO READ THEM ALL IN AND TAKE THE STDEV ACROSS THEM
    fname_end = 'nside%s_sigma%sarcmin_%s.fits' %(nside,sigma,ZBlabel)
    t1 = time.time()
    
    if nside>2048 and Num_Noise==100:
        # you won't be able to read in all noise realisations at once (Memory Error),
        # so read in the 0-25, 25-50, 50-75, 75-100 NoiseKappaStd maps you've already made
        # and add them in quadrature to get the NoiseKappaStd for 0-100
        print("Reading in the pre-made 0-25, 25-50, 50-75 & 75-100 NoiseKappaStd maps.")
        kappa_rng_1 = hp.read_map('MassMaps/K%s_NoiseKappaStd0-25_%s' %(NorS,fname_end))
        kappa_rng_2 = hp.read_map('MassMaps/K%s_NoiseKappaStd25-50_%s' %(NorS,fname_end))
        kappa_rng_3 = hp.read_map('MassMaps/K%s_NoiseKappaStd50-75_%s' %(NorS,fname_end))
        kappa_rng_4 = hp.read_map('MassMaps/K%s_NoiseKappaStd75-100_%s' %(NorS,fname_end))
        
        kappa_rng_std = np.sqrt( kappa_rng_1**2 + kappa_rng_2**2 + kappa_rng_3**2 + kappa_rng_4**2 )

    else:
        kappa_rng_array = np.zeros([ Num_Noise, int(12*nside**2) ])
        k=0
        for i in range( Noise_start, Noise_end ):
            print("Reading in Noise %s of %s" %(i,Noise_end))
            kappa_rng_array[k,:] = hp.read_map('MassMaps/K%s_NoiseKappa%s_%s' %(NorS,i,fname_end))
            k+=1
        kappa_rng_std = np.std( kappa_rng_array, axis=0 )
        #del kappa_rng_array 

    # Write the overall NoiseKappaStd map to file
    hp.write_map('MassMaps/K%s_NoiseKappaStd%s-%s_%s' %(NorS,Noise_start,Noise_end,fname_end),
                 kappa_rng_std, overwrite=True)
    t2 = time.time()
    print( "Reading in %s noise maps took %.1f s" %(Num_Noise,(t2-t1)) )
    
elif Run_Or_Read_Noise == "Read":
    # USE THIS IF THE STDEV ACROSS THE NOISE MAPS IS ALREADY MADE AND JUST NEEDS READING IN.
    kappa_rng_std = hp.read_map('MassMaps/K%s_NoiseKappaStd%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel))
        
else:
    print( "Not running or reading any noise realisations, will read pre-saved SNR map." )


if Read_SNR:
    SNR = hp.read_map('MassMaps/K%s_SNR%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel))
else:
    SNR = kappa_hp / kappa_rng_std
    SNR[ mask ] = hp.UNSEEN
    # Fix holes and save SNR map
    SNR = Fix_Holes(SNR)
    hp.write_map('MassMaps/K%s_SNR%s-%s_nside%s_sigma%sarcmin_%s.fits' %(NorS,Noise_start,Noise_end,nside,sigma,ZBlabel),
             SNR, overwrite=True)


# Add full Moon to blank space in Northern field, as a sense of scale?
Add_Moon = False
if Add_Moon:
    ra_Moon = 150.
    dec_Moon = 0.
    scale = 5 # scale up the radius of the Moon by this factor (and its area by this factor squared)
    print("Adding full Moon to northern field at position RA,Dec %s,%s" %(ra_Moon,dec_Moon))
    radius_Moon = np.radians( 0.5*(31./60.) ) # angular radius of Moon (diameter~31 arcmin in radians)
    SNR_Moon = np.copy(SNR)
    vec = hp.ang2vec( np.radians(dec_Moon)+np.pi/2., np.radians(ra_Moon) ) # convert coordinates to 3D vector
    pix_Moon = hp.query_disc(nside, vec, radius=scale*radius_Moon) #, inclusive=True, fact=4)
    # Sanity checks that the returned pixels are centred on (ra_Moon, dec_Moon):
    tmp_ra_Moon, tmp_dec_Moon = hp.pixelfunc.pix2ang(nside, pix_Moon,lonlat=True)
    SNR_Moon[ pix_Moon ] = 100 # arbitrary high value to stand out.
    



# PRESS RELEASE PLOTS
# NORTH
ras_N = np.arange(140,220,10)
Plot_Orthview((15,2), SNR, [182.,0.], -5, 5,
              'MassMaps/orthview_KN_SNR%s-%s_nside%s_sigma%sarcmin_%s.png' %(Noise_start,Noise_end, nside,sigma,ZBlabel),
              [-0.84,0.82], [-0.1,0.1], ras_N, -10, 10)

# SOUTH
ras_S = np.arange(-20,70,10)
# wider pan plot (for when you're plotting ra&dec labels) [-0.65,0.65], [-0.3,0.25]
#Plot_Orthview((12,4), SNR, [10.,-30.], -5, 5,
#              'MassMaps/orthview_KS_SNR%s-%s_nside%s_sigma%sarcmin_%s.png' %(Noise_start,Noise_end, nside,sigma,ZBlabel),
#              [-0.65,0.6], [-0.22,0.08], ras_S, -40, -20) 

