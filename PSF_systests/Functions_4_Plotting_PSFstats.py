# 15/05/2020, B. M. Giblin, Edinburgh
# Tidying up codes that plot rho/Paulin-Henriksson stats
# by having some of the functions in here.

import numpy as np
from scipy.stats import binned_statistic_2d
from astropy.io import fits
import time

def interpolate2D(X, Y, grid): #(It's linear)                                                                    
    # This function does simple 2D linear interpolation to find values of 'grid' at positions (X,Y)
    # This is faster over, eg, scipy and numpy functions, which have to be performed in a loop.               
    Xi = X.astype(np.int)
    Yi = Y.astype(np.int) # these round down to integer                                                       

    VAL_XYlo = grid[Yi, Xi] + (X - Xi)*( grid[Yi, Xi+1] - grid[Yi, Xi] )
    VAL_XYhi = grid[Yi+1,Xi] + (X - Xi)*( grid[Yi+1,Xi+1] - grid[Yi+1, Xi] )
    VAL_XY = VAL_XYlo + (Y - Yi)*( VAL_XYhi - VAL_XYlo )
    return VAL_XY


def Select_Patch(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the elements in Q corresponding to INSIDE the (ra,dec) range
    idx_ra = np.where(np.logical_and(ra<rhi, ra>rlo))[0]
    idx_dec = np.where(np.logical_and(dec<dhi, dec>dlo))[0]
    idx = np.intersect1d(idx_ra, idx_dec)
    return Q[idx]


def MeanQ_VS_XY(Q, w, X,Y,num_XY_bins):
    # we want the weighted mean of Q.
    # Calculate the sum 2D binned value of Q*w and w, and then divide
    sumQw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q*w, statistic='sum', bins=num_XY_bins)
    sum_w_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, w, statistic='sum', bins=num_XY_bins)
    AvQ_grid=sumQw_grid/sum_w_grid
    return AvQ_grid,sum_w_grid,yedges,xedges

def Bootstrap_Error(nboot, samples, weights):
    N = len(samples)
    bt_samples = np.zeros(nboot)                            # Will store mean of nboot resamples              
    for i in range(nboot):
        idx = np.random.randint(0,N,N)                      # Picks N random indicies with replacement        
        bt_samples[i] = np.sum( weights[idx]*samples[idx] ) /np.sum (weights[idx])
    return np.std( bt_samples )



def Read_GalData(NorS):
    # Use the Master catalogue:                                                                         
    data_DIR = '/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/'
    f = fits.open('%s/K1000_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.cat' %(data_DIR,NorS))
    RA = f[1].data['ALPHA_J2000']
    Dec= f[1].data['DELTA_J2000']
    Z  = f[1].data['Z_B']

    #Calculated the integrals with wolfram                                                                      
    # - turns out we need a factor of 4 to make T_g and T_PSF equivalent                                        
    #https://docs.google.com/document/d/1kylhRNInqzofQtZDiLrH-GH5RwhEBB04OZ0HW8pRZOo/edit?usp=sharing           

    T_gal = 4*f[1].data['bias_corrected_scalelength_pixels']**2
    T_PSF = f[1].data['PSF_Q11'] + f[1].data['PSF_Q22'] # The PSF size at the location of the galaxy            
    weight= f[1].data['recal_weight_A']*f[1].data['Flag_SOM_Fid_A'] #Include SOM Flag in the weight             

    return RA, Dec, T_gal, T_PSF, weight, Z


def Interp_deltaT(RA_g, Dec_g, RA_p, Dec_p, delta_TPSF):
    # This function interpolates delta_TPSF to the galaxy positions
    
    # To grid up the survey, decide on a suitable angular size a pixel should be:                                     
    # Use the dimensions of the PSF data (_p), not gal data (_g),                                                     
    # as the PSF data spans wider (RA,Dec)                                                                            
    ang_pxl = 5. / 60. # 5 arcmin in degrees                                                                          
    nbins_x = int( ( RA_p.max()-RA_p.min() ) / ang_pxl )
    nbins_y = int( ( Dec_p.max()-Dec_p.min() ) / ang_pxl )
    # pixel coordinates of the PSF objects                                                                            
    X_p = nbins_x * (RA_p - RA_p.min()) / (RA_p.max()-RA_p.min())
    Y_p = nbins_y * (Dec_p - Dec_p.min()) / (Dec_p.max()-Dec_p.min())

    # Turn delta_TPSF into a grid:                                                                                    
    delta_TPSF_grid,count_grid, _,_ = MeanQ_VS_XY(delta_TPSF, np.ones_like(delta_TPSF),
                                                  X_p,Y_p, [nbins_y,nbins_x])
    delta_TPSF_grid = np.nan_to_num( delta_TPSF_grid, nan=0. ) # Lots of nans due to 0/0                              
    # Need to append TPSF_grid with final row and column to avoid interp error                                        
    delta_TPSF_grid = np.c_[ delta_TPSF_grid, delta_TPSF_grid[:,-1] ]
    delta_TPSF_grid = np.r_[ delta_TPSF_grid, [delta_TPSF_grid[-1,:]] ]

    # pixel coordinates of galaxies                                                                                   
    X_g = nbins_x * (RA_g - RA_p.min()) / (RA_p.max()-RA_p.min())
    Y_g = nbins_y * (Dec_g - Dec_p.min()) / (Dec_p.max()-Dec_p.min())

    # Finally get delta_T_PSF at the position of the galaxies                                                         
    delta_TPSF_g = interpolate2D(X_g, Y_g, delta_TPSF_grid)
    return delta_TPSF_g



def Calc_Important_Tquantities(LFver,zbounds, nboot):
    # This functions calculates:
    # < deltaT_PSF / T_gal > & SIGMA[ deltaT_PSF / T_gal ]
    # < T_gal^-2 > & SIGMA[ T_gal^-2  ] 

    # Read in N and S PSF data catalogues
    # '_p' means at position of objects used for PSF modelling
    bgdir='/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/PSF_systests'
    RA_Np, Dec_Np, _,_,_,_, TPSF_Np, delta_TPSF_Np, _,_ = np.load('%s/LFver%s/Catalogues/PSF_Data_N.npy'%(bgdir,LFver[0])).transpose()
    RA_Sp, Dec_Sp, _,_,_,_, TPSF_Sp, delta_TPSF_Sp, _,_ = np.load('%s/LFver%s/Catalogues/PSF_Data_S.npy'%(bgdir,LFver[0])).transpose()

    # Read in N and S galaxy data catalogues
    RA_Ng, Dec_Ng, T_Ng, T_PSF_Ng, weight_Ng, z_Ng = Read_GalData('N')
    RA_Sg, Dec_Sg, T_Sg, T_PSF_Sg, weight_Sg, z_Sg = Read_GalData('S')

    # Now scroll through the z-bins, calculating the T-quantities in each.
    deltaT_ratio = np.zeros([ 2, len(zbounds)-1 ]) # 1st-row = mean, 2nd-row = error
    Tg_invsq = np.zeros_like( deltaT_ratio )       # Same
    for i in range(len(zbounds)-1):
        t1 = time.time()
        print('On tomo bin %s-%s' %(zbounds[i],zbounds[i+1]))
        # Redshift cut:
        idx_N = ( (z_Ng>zbounds[i]) & (z_Ng<zbounds[i+1]) )
        idx_S = ( (z_Sg>zbounds[i]) & (z_Sg<zbounds[i+1]) )
        delta_TPSF_Ng = Interp_deltaT(RA_Ng[idx_N], Dec_Ng[idx_N], RA_Np, Dec_Np, delta_TPSF_Np)
        delta_TPSF_Sg = Interp_deltaT(RA_Sg[idx_S], Dec_Sg[idx_S], RA_Sp, Dec_Sp, delta_TPSF_Sp)

        # 1. < deltaT_PSF / T_gal > & SIGMA[ deltaT_PSF / T_gal ]
        # NOTE: T_g is on the DENOMINATOR here, so don't think weights shoulf be weight_g !
        weight_g = np.append(weight_Ng[idx_N], weight_Sg[idx_S]) 
        dT_p = np.append(delta_TPSF_Ng, delta_TPSF_Sg)
        T_g = np.append(T_Ng[idx_N], T_Sg[idx_S])
        deltaT_ratio[0,i] = np.average( dT_p/T_g, weights=weight_g )
        # This is a weighted mean, so lets use a bootstrap estimate for the error on the mean
        print("Bootstrapping deltaT-ratio with nboot=%s" %nboot)
        deltaT_ratio[1,i] = Bootstrap_Error(nboot, dT_p/T_g, weights=weight_g)

        # Query, what should the weights be here?
        # 2. < T_gal^-2 > & SIGMA[ T_gal^-2  ] 
        Tg_invsq[0,i] = np.average( 1./T_g**2., weights=weight_g )
        print("Bootstrapping Tgal_invsq with nboot=%s" %nboot)
        Tg_invsq[1,i] = Bootstrap_Error(nboot, 1/T_g**2., weights=weight_g )

        t2 = time.time()
        print('For tomo bin %s-%s, got the following T-quantities (took %0.f s)' %(zbounds[i],zbounds[i+1],(t2-t1)) )
        print ('%8.3e,%8.3e,%8.3e,%8.3e'%(deltaT_ratio[0,i], deltaT_ratio[1,i], Tg_invsq[0,i], Tg_invsq[1,i]))

    return deltaT_ratio, Tg_invsq

