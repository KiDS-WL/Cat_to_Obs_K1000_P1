# 15/05/2020, B. M. Giblin, Edinburgh
# Tidying up codes that plot rho/Paulin-Henriksson stats
# by having some of the functions in here.

import numpy as np
from scipy.stats import binned_statistic_2d
from astropy.io import fits
import time
import glob

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
    
    num_zbins = len(zbounds)-1
    num_zbins_tot = np.sum( range(num_zbins+1) ) # Includes cross bins
    
    # Read in N and S PSF data catalogues
    # '_p' means at position of objects used for PSF modelling
    bgdir='/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/PSF_systests'
    RA_Np, Dec_Np, _,_,_,_, TPSF_Np, delta_TPSF_Np, _,_ = np.load('%s/LFver%s/Catalogues/PSF_Data_N.npy'%(bgdir,LFver)).transpose()
    RA_Sp, Dec_Sp, _,_,_,_, TPSF_Sp, delta_TPSF_Sp, _,_ = np.load('%s/LFver%s/Catalogues/PSF_Data_S.npy'%(bgdir,LFver)).transpose()

    # Read in N and S galaxy data catalogues
    RA_Ng, Dec_Ng, T_Ng, T_PSF_Ng, weight_Ng, z_Ng = Read_GalData('N')
    RA_Sg, Dec_Sg, T_Sg, T_PSF_Sg, weight_Sg, z_Sg = Read_GalData('S')

    # The fact that the S RA data crosses zero causes issues interpolating onto a grid.
    # So shift rescale all S RA's in [300,360] to be negative, making the field continuous.
    RA_Sg[ ((RA_Sg<360) & (RA_Sg>300)) ] += -360.
    RA_Sp[ ((RA_Sp<360) & (RA_Sp>300)) ] += -360.
    
    # Now scroll through the z-bins, calculating the T-quantities in each.
    deltaT_ratio = np.zeros([ 2, num_zbins ])      # 1st-row = mean, 2nd-row = error
    Tg_invsq = np.zeros_like( deltaT_ratio )       # Same
    for i in range(num_zbins):
        t1 = time.time()
        print('On tomo bin %s-%s' %(zbounds[i],zbounds[i+1]))
        # Redshift cut:
        idx_N = ( (z_Ng>zbounds[i]) & (z_Ng<zbounds[i+1]) )
        idx_S = ( (z_Sg>zbounds[i]) & (z_Sg<zbounds[i+1]) )
        delta_TPSF_Ng = Interp_deltaT(RA_Ng[idx_N], Dec_Ng[idx_N], RA_Np, Dec_Np, delta_TPSF_Np)
        delta_TPSF_Sg = Interp_deltaT(RA_Sg[idx_S], Dec_Sg[idx_S], RA_Sp, Dec_Sp, delta_TPSF_Sp)

        # 1. < deltaT_PSF / T_gal > & SIGMA[ deltaT_PSF / T_gal ]
        # NOTE: T_g is on the DENOMINATOR here, so don't think weights shoulf be weight_g !
        T_g = np.append(T_Ng[idx_N], T_Sg[idx_S])
        dT_p = np.append(delta_TPSF_Ng, delta_TPSF_Sg)
        weight_g = np.append(weight_Ng[idx_N], weight_Sg[idx_S])
        weight_q1 = weight_g #(T_g**2 * weight_g ) / dT_p
        #weight_q1 = np.nan_to_num( weight_q1, nan=0., posinf=0. ) # get rid of nan/inf weights from where dT_p=0
        # ^ that weight computed with error propagation:
        # if y=a/x, sigma_y^2=(-a/x^2)sigma_x^2, sigma_(x/y)^2=1/weight_(x/y)
        # then x=T_g, a=dT_p, y=dT_p/T_g
        
        deltaT_ratio[0,i] = np.average( dT_p/T_g, weights=weight_q1 )
        # This is a weighted mean, so lets use a bootstrap estimate for the error on the mean
        print("Bootstrapping deltaT-ratio with nboot=%s" %nboot)
        deltaT_ratio[1,i] = Bootstrap_Error(nboot, dT_p/T_g, weights=weight_q1)

        # 2. < T_gal^-2 > & SIGMA[ T_gal^-2  ]
        weight_q2 = weight_g  #weight_g #T_g**2 * weight_g
        # ^ computed with error prop: y=1/x, sigma_y^2=(dy/dx)sigma_x^2, sigma_x^2=1/weight_x
        # with x=T_g, weight_x=weight_g
        Tg_inv = np.average( 1./T_g, weights=weight_q2 )
        Tg_invsq[0,i] = Tg_inv**2
        print("Bootstrapping Tgal_invsq with nboot=%s" %nboot)
        Tg_inverr = Bootstrap_Error(nboot, 1/T_g, weights=weight_q2 )
        # Need to convert ^this error on 1/T_g to an error on 1/T_g^2
        # do error propagation again:
        # z=y^2, sigma_z=2y*sigma_y^2
        Tg_invsq[1,i] = np.sqrt(2 * Tg_inv) * Tg_inverr

        t2 = time.time()
        print('For tomo bin %s-%s, got the following T-quantities (took %0.f s)' %(zbounds[i],zbounds[i+1],(t2-t1)) )
        print ('%8.3e,%8.3e,%8.3e,%8.3e'%(deltaT_ratio[0,i], deltaT_ratio[1,i], Tg_invsq[0,i], Tg_invsq[1,i]))

    # For the cross-bins, for now just average the T-quantities in the individual bins
    deltaT_ratio_tot = np.zeros([ 2, num_zbins_tot ])      # 1st-row = mean, 2nd-row = error
    Tg_invsq_tot = np.zeros_like( deltaT_ratio_tot )           # Same
    k=0
    for i in range(num_zbins):
        for j in range(num_zbins):
            if j>= i:
                deltaT_ratio_tot[0,k] = (deltaT_ratio[0,i]+deltaT_ratio[0,j])/2
                deltaT_ratio_tot[1,k] = np.sqrt( deltaT_ratio[1,i]**2 + deltaT_ratio[1,j]**2 )
                Tg_invsq_tot[0,k] = (Tg_invsq[0,i] + Tg_invsq[0,j])/2
                Tg_invsq_tot[1,k] = np.sqrt(Tg_invsq[0,i]**2 + Tg_invsq[0,j]**2)
                k+=1
                
    return deltaT_ratio, Tg_invsq, deltaT_ratio_tot, Tg_invsq_tot



def Read_rho_Or_PH(LFver, keyword, ThBins, Res):
    # This either reads in and returns the rho stats or the Paulin-Henriksson terms
    # depending on what keyword is set to.
    # Note the unfortunate variable naming here 'ph' (Paulin-Henriksson)
    # even in the case of rho stats being read in.
    
    if keyword == 'rho':
        print("Reading in the rho stats")
        num_stat = 5
        DIR=['rho1','rho2','rho3','rho4','rho5']
        stat='rho'
    else:
        print("Reading in the Paulin-Henriksson terms")
        num_stat = 3
        DIR=['PHterms','PHterms','PHterms','PHterms','PHterms']
        stat='ph'

    NFiles = []
    numN = []
    SFiles = []
    numS = []
    Plot_Labels = []
    for lfv in range(len(LFver)):
        NFiles.append( glob.glob('LFver%s/%s/%s1_KiDS_N_*of%sx%s.dat'%(LFver[lfv],DIR[0],stat,Res,Res)) )
        numN.append(len( NFiles[lfv] ))
        SFiles.append( glob.glob('LFver%s/%s/%s1_KiDS_S_*of%sx%s.dat'%(LFver[lfv],DIR[0],stat,Res,Res)) )
        numS.append(len( SFiles[lfv] ))
        if LFver[lfv] == "309b":
            Plot_Labels.append(" 3:1")
        elif LFver[lfv] == "319":
            Plot_Labels.append(LFver[lfv] + " 3:1")
        elif LFver[lfv] == "319b":
            Plot_Labels.append(" 4:1")
        elif LFver[lfv] == "319c":
            Plot_Labels.append(" 3:2")
        elif LFver[lfv] == "319d":
            Plot_Labels.append(" 5:1")
        elif LFver[lfv] == "321":
            Plot_Labels.append("New 4:1")

        
    php_mean = np.zeros([len(LFver),num_stat,ThBins])
    php_err = np.zeros_like(php_mean)
    # Fix theta to the values saved for the very first statistic:
    theta = np.loadtxt('LFver%s/%s/%s%s_KiDS_N.dat'%(LFver[0],DIR[0],stat,1), usecols=(0,), unpack=True)
    
    for lfv in range(len(LFver)):

        php_split = np.zeros([ numN[lfv]+numS[lfv], 5, ThBins ])

        for i in range(num_stat):
            try:
                tmp_theta, phpN = np.loadtxt('LFver%s/%s/%s%s_KiDS_N.dat'%(LFver[lfv],DIR[i],stat,i+1), usecols=(0,1), unpack=True)
                # If the above exists, try to read in the weight (only saved this for LFver321 of the rho stats)
                try:
                    weightN = np.loadtxt('LFver%s/%s/%s%s_KiDS_N.dat'%(LFver[lfv],DIR[i],stat,i+1), usecols=(3,), unpack=True)
                except IndexError:
                    weightN = 1.
            except IOError:
                weightN = np.ones(len(ThBins))
                phpN = np.zeros(len(ThBins))

            # Interp to the fixed_theta scale
            phpN = np.interp( theta, tmp_theta, phpN )
            # keep the weights as they are - don't interp.
                
            try:
                tmp_theta, phpS = np.loadtxt('LFver%s/%s/%s%s_KiDS_S.dat'%(LFver[lfv],DIR[i],stat,i+1), usecols=(0,1), unpack=True)
                # If the above exists, try to read in the weight (only saved this for LFver321 of the rho stats)
                try:
                    weightS = np.loadtxt('LFver%s/%s/%s%s_KiDS_S.dat'%(LFver[lfv],DIR[i],stat,i+1), usecols=(3,), unpack=True)
                except IndexError:
                    weightS = 1.
            except IOError:
                weightS = np.ones(len(ThBins))
                phpS = np.zeros(len(ThBins))

            # Interp to the fixed_theta scale
            phpS = np.interp( theta, tmp_theta, phpS )
            # keep the weights as they are - don't interp. 
                
            php_mean[lfv,i,:] = (weightN*phpN + weightS*phpS) / (weightN+weightS)
            for j in range(numN[lfv]):
                tmp_theta, tmp_php_splitN = np.loadtxt(NFiles[lfv][j], usecols=(0,1), unpack=True)
                php_split[j,i,:] = np.interp( theta, tmp_theta, tmp_php_splitN )
                
            for j in range(numS[lfv]):
                tmp_theta, tmp_php_splitS = np.loadtxt(SFiles[lfv][j], usecols=(0,1), unpack=True)
                php_split[numN[lfv]+j,i,:] = np.interp( theta, tmp_theta, tmp_php_splitS )

            php_err[lfv,i,:] = np.sqrt( np.diag( np.cov(php_split[:,i,:], rowvar = False) ) / (numN[lfv]+numS[lfv]) )

    return Plot_Labels, theta, php_mean, php_err




def Read_alpha_per_bin(LFver):
    # This functions reads in the alpha (PSF leakage) values for the 5 tomo bins
    # and avg's to estimate values for the cross-bins,

    num_zbins = 5
    num_zbins_tot = np.sum( range(num_zbins+1) )    # Number source bins including cross-bins
    alpha = np.zeros([ len(LFver), num_zbins_tot ])
    for lfv in range(len(LFver)):

        if LFver[lfv] == "321":
            tmp = "glab_%s"%LFver[lfv]
        elif LFver[lfv] == "309c":
            tmp = "svn_%s"%LFver[lfv]
        else:
            print("Currently only have saved alpha values for 2 LF versions: 321 and 309c. EXITING")
            sys.exit()

        tmp_a1, tmp_a2 = np.loadtxt('KAll.autocal.BlindA.alpha_VS_ZB.ZBcut0.1-1.2_LF_%s_2Dbins.dat'%tmp,
                                    usecols=(1,3), unpack=True)
        tmp_a = (tmp_a1 + tmp_a2) / 2.    # Just taking the avg of the alpha per ellipticity component                                 
                                          # Also calculate the alphas in the cross-bins as the combination of values in the auto-bins  
        k=0
        for i in range(num_zbins):
            for j in range(num_zbins):
                if j>= i:
                    alpha[lfv,k] = np.sqrt( tmp_a[i]*tmp_a[j] )
                    #print("%s : %s %s : %s" %(k, tmp_a1[i], tmp_a1[j], alpha[lfv,k]) )                                    
                    k+=1

    return alpha


def Read_In_Theory_Vector(hi_lo_fid):
    num_zbins_tot = 15
    # This function returns the theoretical xi+ predictions with either a fiducial, high or low S8 value
    # hi_lo_fid must be one of 'high', 'low' or 'fid'

    # THESE LINES PULL FROM THE OLD THEORY PREDICTIONS (KV450 nofz)                                                                        
    #indir_theory = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/outputs/test_output_S8_%s_test/shear_xi_plus/' %hi_lo_fid                                                                                                                                          
    #theta_theory = np.loadtxt('%s/theta.txt' %indir_theory) * (180./np.pi) * 60.   # Convert long theta array in radians to arcmin        

    indir_theory = '/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/new_outputs/test_output_S8_%s_test/chain/output_test_A/shear_xi_plus_binned/' %hi_lo_fid
    theta_theory = np.loadtxt('%s/theta_bin_1_1.txt' %indir_theory)
    xip_theory_stack = np.zeros( [num_zbins_tot,len(theta_theory)] )
    # ^This will store all auto & cross xi_p for the 5 tomo bins                                                                           
    idx = 0
    for i in range(1,6):
        for j in range(1,6):
            if i >= j:              # Only read in bins 1-1, 2-1, 2-2, 3-1, 3-2,...                                                
                tmp_theta_theory = np.loadtxt('%s/theta_bin_%s_%s.txt' %(indir_theory,i,j))
                tmp_xip_theory = np.loadtxt('%s/bin_%s_%s.txt' %(indir_theory,i,j))
                xip_theory_stack[idx,:] = np.interp( theta_theory, tmp_theta_theory, tmp_xip_theory )
                # pretty sure theta_theory is the same for every tomo bin combination
                # but just in case, we interpolate everything to be sampled at the values for the 1_1 bin.
                idx+=1
    xip_theory_stack = xip_theory_stack.flatten()
    return theta_theory, xip_theory_stack  
