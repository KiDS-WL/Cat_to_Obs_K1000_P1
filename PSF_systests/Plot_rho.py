import numpy as np
import pylab as plt
from astropy.io import fits
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import chi2
from scipy.stats import multivariate_normal as multi_norm
from scipy.stats import binned_statistic_2d
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.stats import ks_2samp	# For KS tests
import time

import sys
import os 
import glob

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14} # 19

plt.rc('font', **font)


# NOTE: THE FOLLOWING VARIABLES ARE USED TO SET REQUIREMENTS ON RHO-STATS FOR
# VISUAL PURPOSES ONLY.
# ------------------------------------------------ #
# Set the plotted requirement on the rho stats:
Requirement = "M18"			# "M18" for Mandelbaum+18: [p_1,2,4<xi_+/SNR, p_2,5<xi_+/(SNR*alpha)
							# "Z18" for Zuntz+18: xi_+ / 10 

Which_Data = "Data"			# "Data" for K1000 or "Theory" for HaloFit
							# goes into calculation of the requirement bands.

# If using the lowest ZB bin, you might want to smooth the data vector 
# so the Requirement band is less jaggardy. If so, use this smoothing scale:
Smooth_Scale = 0.8
# ------------------------------------------------ #

# Define quantities needed for delta_xip
LFver = ["321"] #["309b",  "319b", "319c", "319d", "321"] # "319",

T_ratio = 0.9976 	 # average T_PSF / T_gal Ratio
T_ratio_err = 0.0795
deltaT_ratio= -0.0012    # average deltaT_PSF / T_PSF
deltaT_ratio_err = 0.0263
                         # These numbers are quite accurate, but use 'Calculate_Important_Tquantities()'
                         # below to calc them directly from the catalogues if you want proof.


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
        # we want the weighted mean of Q and also calibrated.  Calculate the sum 2D binned value of Q*w and m*w,
        # and then divide
        sumQw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q*w, statistic='sum', bins=num_XY_bins)
        sum_w_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, w, statistic='sum', bins=num_XY_bins)
        AvQ_grid=sumQw_grid/sum_w_grid
        return AvQ_grid,sum_w_grid,yedges,xedges


def Bootstrap_Error(nboot, samples, weights):
	N = len(samples)
	bt_samples = np.zeros(nboot)		 		# Will store mean of nboot resamples
	for i in range(nboot):
		idx = np.random.randint(0,N,N)			# Picks N random indicies with replacement
		bt_samples[i] = np.sum( weights[idx]*samples[idx] ) /np.sum (weights[idx])
	return np.std( bt_samples )

def Calc_Important_Tquantities(zmin,zmax):
        # Read in N and S catalogues
        
        bgdir='/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/PSF_systests'
        RA_N, Dec_N, e0PSF_N, e1PSF_N, delta_e0PSF_N, delta_e1PSF_N, TPSF_N, delta_TPSF_N, Xpos_N, Ypos_N = np.load('%s/LFver%s/Catalogues/PSF_Data_N.npy'%(bgdir,LFver[0])).transpose()
        RA_S, Dec_S, e0PSF_S, e1PSF_S, delta_e0PSF_S, delta_e1PSF_S, TPSF_S, delta_TPSF_S, Xpos_S, Ypos_S = np.load('%s/LFver%s/Catalogues/PSF_Data_S.npy'%(bgdir,LFver[0])).transpose()

        # Append them together
        RA = np.append(RA_N, RA_S)
        Dec = np.append(Dec_N, Dec_S)
        delta_e0PSF = np.append( delta_e0PSF_N, delta_e0PSF_S )
        delta_e1PSF = np.append( delta_e1PSF_N, delta_e1PSF_S )
        TPSF = np.append( TPSF_N, TPSF_S )
        delta_TPSF = np.append( delta_TPSF_N, delta_TPSF_S )
        
	# Calc important quantities:
	# First < delta T_PSF / T_PSF >
        delta_Tratio = np.mean( delta_TPSF / TPSF )      # Comes out =-0.00128
        # This gives the variance of the PSF residual size ratio
        # delta_Tratio_err = np.std( delta_TPSF / TPSF )   
        # but we want the error on the mean measurement
        delta_Tratio_err = np.std( delta_TPSF / TPSF )/np.sqrt(len(TPSF))   # Comes out = 7.8e-6
        
        # Secondly: < T_PSF / T_gal > and < delta T_PSF / Tgal>
        # Need to interpolate delta_T_PSF to galaxy positions - read in (RA,Dec) of sources
        def Read_GalData(NorS):
                # This read the old nofz version of the catalogue.
                #data_g = np.load('../Calc_1pt_Stats/Catalogues/K%s.BlindA.ra_dec_e1_e2_w_ZB.ZBcutNone.npy'%NorS)
                #RA_g = data_g[:,0]
                #Dec_g = data_g[:,1]
                #Z_g = data_g[:,5]

                # Instead use the Master catalogue:
                data_DIR = '/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/'
                f = fits.open('%s/K1000_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.cat' %(data_DIR,NorS))
                RA_g = f[1].data['ALPHA_J2000']
                Dec_g= f[1].data['DELTA_J2000']
                Z_g  = f[1].data['Z_B']
                
                #From 3.4.1: We estimate the average unconvolved galaxy size Tgal, by taking the lensfit- weighted average of the exponential disk scalelength (squared),
                #as determined from each best-fit galaxy model. We note that this choice results in inconsistent definitions for PSF size and galaxy size, but
                #argue this is an improvement on the alternative approach of setting TPSF/Tgal = 1 

                T_g = f[1].data['bias_corrected_scalelength_pixels']**2
                T_PSF = f[1].data['PSF_Q11'] + f[1].data['PSF_Q22'] # The PSF size at the location of the galaxy
                weight= f[1].data['recal_weight_A']*f[1].data['Flag_SOM_Fid_A'] #Include SOM Flag in the weight 
                
                #idx = ( (Z_g>0.1) & (Z_g<1.2) ) # Redshift cut - need to calculate this for the different tomo bins
                idx = ( (Z_g>zmin) & (Z_g<zmax) ) # Redshift cut - need to calculate this for the different tomo bins
                return RA_g[idx], Dec_g[idx], T_g[idx], T_PSF[idx], weight[idx]
        RA_Ng, Dec_Ng, T_Ngal, T_NPSF, weight = Read_GalData('N')
	#RA_Sg, Dec_Sg = Read_GalData('S') # Haven't saved the S catalogue apparently
        
        # To grid up the survey, decide on a suitable angular size a pixel should be:
        # Use the dimensions of the PSF DATA, not gal data, as PSF spans wider (RA,Dec)
        ang_pxl = 5. / 60. # 5 arcmin in degrees
        nbins_x = int( ( RA_N.max()-RA_N.min() ) / ang_pxl )
        nbins_y = int( ( Dec_N.max()-Dec_N.min() ) / ang_pxl )
	# pixel coordinates of the PSF objects
        X_N = nbins_x * (RA_N - RA_N.min()) / (RA_N.max()-RA_N.min()) 
        Y_N = nbins_y * (Dec_N - Dec_N.min()) / (Dec_N.max()-Dec_N.min()) 
        # Ammend the value at MAX(X,Y) which will be at edge of map and cause interp error.
        #X_N[np.argmax(X_N)] += -0.1
        #Y_N[np.argmax(Y_N)] += -0.1
        # Turn TPSF into a grid:
        delta_TPSF_grid,count_grid, _,_ = MeanQ_VS_XY(delta_TPSF_N, np.ones_like(delta_TPSF_N), X_N,Y_N, [nbins_y,nbins_x])
        delta_TPSF_grid = np.nan_to_num( delta_TPSF_grid, nan=0. ) # Lots of nans due to 0/0
        # Need to append TPSF_grid with final row and column to avoid interp error
        delta_TPSF_grid = np.c_[ delta_TPSF_grid, delta_TPSF_grid[:,-1] ] 
        delta_TPSF_grid = np.r_[ delta_TPSF_grid, [delta_TPSF_grid[-1,:]] ] 

        # pixel coordinates of galaxies
        X_Ng = nbins_x * (RA_Ng - RA_N.min()) / (RA_N.max()-RA_N.min()) 
        Y_Ng = nbins_y * (Dec_Ng - Dec_N.min()) / (Dec_N.max()-Dec_N.min()) 
	
        # Finally get delta_T_PSF at the position of the galaxies
        delta_TPSF_atgal = interpolate2D(X_Ng, Y_Ng, delta_TPSF_grid)

        T_ratio = np.average( T_NPSF/ T_Ngal, weights=weight)    #
        # This is a weighted mean, so lets use a bootstrap estimate for the error on the mean
        #T_ratio_err = np.std( T_NPSF / T_Ngal ) #
        nboot=30
        T_ratio_err=Bootstrap_Error(nboot, T_NPSF/T_Ngal, weight)
        
        delta_TPSF_over_Tgal = np.average( delta_TPSF_atgal/ T_Ngal, weights=weight)    # 
        delta_TPSF_over_Tgal_err = Bootstrap_Error(nboot, delta_TPSF_atgal/T_Ngal, weight)
        print ('%8.3e,%8.3e,%8.3e,%8.3e'%(T_ratio, T_ratio_err, delta_TPSF_over_Tgal, delta_TPSF_over_Tgal_err))
        return delta_Tratio, delta_Tratio_err, T_ratio, T_ratio_err, delta_TPSF_over_Tgal, delta_TPSF_over_Tgal_err
delta_Tratio, delta_Tratio_err, T_ratio, T_ratio_err, delta_TPSF_over_Tgal, delta_TPSF_over_Tgal_err= Calc_Important_Tquantities(0.1,0.3)
# TO DO 1 (BG) - carry through tomography-dependent ratio values through to the final result
# Bins 2,3,4 look pretty similar, but 1 and 5 are quite different
# TO DO 2 (BG) - add in the South - the PSF is quite different in the South - North on its own isn't fully representative
# TO DO 3 (CH) - how much of a mismatch is Q11+Q22 vs exponetial scale length - do we need a normalisation constant?

# Read in the alpha values for each shear component and tomo-bin
Use_alpha_per_bin = True                        # If True, use an alpha per bin
num_zbins   = 5                                 # Number of source bins
num_zbins_tot = np.sum( range(num_zbins+1) )    # Number source bins including cross-bins
alpha = np.zeros([ len(LFver), num_zbins_tot ])
if Use_alpha_per_bin:
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

else:
	alpha += 0.03     # worse case scenario PSF leakage is 0.03;

ThBins = 9
Res = 7


# Find number of rho_pm contributing to cov
NFiles = []
numN = []
SFiles = []
numS = []
Plot_Labels = []
for lfv in range(len(LFver)):
	NFiles.append( glob.glob('LFver%s/rho1/rho1_KiDS_N_*of%sx%s.dat'%(LFver[lfv],Res,Res)) )
	numN.append(len( NFiles[lfv] ))
	SFiles.append( glob.glob('LFver%s/rho1/rho1_KiDS_S_*of%sx%s.dat'%(LFver[lfv],Res,Res)) )
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

# Read in avg rho's and those used to calc cov
rhop_mean = np.zeros([len(LFver),5,ThBins])
rhop_err = np.zeros_like(rhop_mean)


for lfv in range(len(LFver)):

	rhop_split = np.zeros([ numN[lfv]+numS[lfv], 5, ThBins ]) 
	#meanfactor = 2. # Divide by this, unless no N or S Field, in which case set to 1.
	for i in range(5):
		try:
			theta, rhopN = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_N.dat'%(LFver[lfv],i+1,i+1), usecols=(0,1), unpack=True)
			# If the above exists, try to read in the weight (only saved this for LFver321)
			try:
				weightN = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_N.dat'%(LFver[lfv],i+1,i+1), usecols=(3,), unpack=True)
			except IndexError:
				weightN = 1.
		except IOError:
			weightN = 1.
			rhopN = 0.
		try:
			theta, rhopS = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_S.dat'%(LFver[lfv],i+1,i+1), usecols=(0,1), unpack=True)
			# If the above exists, try to read in the weight (only saved this for LFver321)
			try:
				weightS = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_S.dat'%(LFver[lfv],i+1,i+1), usecols=(3,), unpack=True)
			except IndexError:
				weightS = 1.
		except IOError:
			weightS = 1.
			rhopS = 0.
	
		# Weighted average of rho-North & rho-South
		rhop_mean[lfv,i,:] = (weightN*rhopN + weightS*rhopS) / (weightN+weightS)
		for j in range(numN[lfv]):
			rhop_split[j,i,:] = np.loadtxt(NFiles[lfv][j], usecols=(1,), unpack=True)
		for j in range(numS[lfv]):
			rhop_split[numN[lfv]+j,i,:] = np.loadtxt(SFiles[lfv][j], usecols=(1,), unpack=True)

		rhop_err[lfv,i,:] = np.sqrt( np.diag( np.cov(rhop_split[:,i,:], rowvar = False) ) / (numN[lfv]+numS[lfv]) ) 


# To get a rough idea of size of rho stats, read in the xi+- of some data to overplot 
data_dir = 'xi_pm_Vectors_4_Overplotting/'  # '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/xi_pm/'
data_ZBlabel = '0.7-0.9'					# ZBcut on the vector you overplot... '0.7-0.9', '0.1-0.3' '0.1-2.0'
theta_data, xip_data = np.loadtxt('%s/KAll.BlindA.xi_pm.ZBcut%s.dat' %(data_dir, data_ZBlabel), usecols=(0,1), unpack=True)

# Read in a covariance 
Linc_Rescale = 600. / 878.83	# Linc says I should rescale his cov by approx. this factor
								# to get the effective area right. This probs isn't 100% accurate.
Cov_inDIR = './Lincs_CovMat'
            # eday address: '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/Lincs_CovMat/'
Cov_Mat_uc_Survey = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[81:90, 81:90] * Linc_Rescale * 0.16		
																	# [0:9, 0:9] This extracts xi+ Cov in lowest bin
																	# [81:90, 81:90] This pulls out the middle bin: 3-3

def Read_In_Theory_Vector(hi_lo_fid):
	# hi_lo_fid must be one of 'high', 'low' or 'fid'
	#indir_theory = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/outputs/test_output_S8_%s_test/shear_xi_plus/' %hi_lo_fid
	#theta_theory = np.loadtxt('%s/theta.txt' %indir_theory) * (180./np.pi) * 60.	# Convert long theta array in radians to arcmin

	indir_theory = '/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/new_outputs/test_output_S8_%s_test/chain/output_test_A/shear_xi_plus_binned/' %hi_lo_fid
        # '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/new_outputs/test_output_S8_%s_test/chain/output_test_A/shear_xi_plus_binned/' %hi_lo_fid
	xip_theory_stack = np.zeros( [num_zbins_tot,len(theta_data)] )									# Will store all auto & cross xi_p for the 5 tomo bins
	idx = 0
	for i in range(1,6):
		for j in range(1,6):
			if i >= j:		# Only read in bins 1-1, 2-1, 2-2, 3-1, 3-2,...
				tmp_theta_theory = np.loadtxt('%s/theta_bin_%s_%s.txt' %(indir_theory,i,j)) 
				tmp_xip_theory = np.loadtxt('%s/bin_%s_%s.txt' %(indir_theory,i,j))
				xip_theory_stack[idx,:] = np.interp( theta_data, tmp_theta_theory, tmp_xip_theory )		# sample at the theta values used for PSF modelling.
				idx+=1
	xip_theory_stack = xip_theory_stack.flatten()
	return xip_theory_stack									# [126:135, 126:135] for xi+ in highest bin


def Calc_delta_xip_J16(alphas, T_ratio, rho, rho_err): 		# !!! Jarvis (2016) expression for PSF systematic
	xip_theory = Read_In_Theory_Vector('fid')	
	xip_theory = np.reshape(xip_theory, (num_zbins_tot,len(theta))) # Reshape to be [num_zbins,ntheta]

	# Calculate the additive shear bias for each lensfit version and in each tomographic bin (ie- value of alpha)
	delta_xip = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
	err_delta_xip = np.zeros_like( delta_xip )
	for lfv in range(len(LFver)):
		for j in range(num_zbins_tot):
			delta_xip[lfv,j,:] = 2*xip_theory[j,:]*T_ratio*deltaT_ratio + T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alphas[lfv,j]*(rho[lfv,1,:]+rho[lfv,4,:])
			err_term1 = T_ratio**2*(rho_err[lfv,0,:]**2+rho_err[lfv,2,:]**2+rho_err[lfv,3,:]**2)
			err_term2 = T_ratio*alphas[lfv,j]*(rho_err[lfv,1,:]**2+rho_err[lfv,4,:]**2)
			err_delta_xip[lfv,j,:] = ( err_term1 - err_term2 )**0.5
	return delta_xip, err_delta_xip


def Calc_delta_xip_H20(T_ratio, rho, rho_err): 		# !!! Heymans' (2020) derivation for PSF systematic
	xip_theory = Read_In_Theory_Vector('fid')	
	xip_theory = np.reshape(xip_theory, (num_zbins_tot,len(theta))) # Reshape to be [num_zbins,ntheta]

	# Calculate the additive shear bias for each lensfit version and in each tomographic bin (ie- value of alpha)
	delta_xip_total = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
	delta_xip_terms = np.zeros([ len(LFver), num_zbins_tot, ThBins, 4 ]) # Store the separate ingredients of the total delta_xip

	err_delta_xip = np.zeros_like( delta_xip_total )
	for lfv in range(len(LFver)):
		for j in range(num_zbins_tot):
			delta_xip_terms[lfv,j,:,0] = 2*xip_theory[j,:]*delta_TPSF_over_Tgal
                        #TO DO 4 (BG) - we need to update the rho calculation in so the prefactor is 1/T_gal^2 rather than T_ratio**2
                        #The current set-up ignores potential cross correlation between the residuals and PSF Size
                        #" we choose to keep all terms that may couple within the correlation function"
                        delta_xip_terms[lfv,j,:,1] = T_ratio**2 *(rho[lfv,0,:])
			delta_xip_terms[lfv,j,:,2] = T_ratio**2 *(rho[lfv,2,:])
			delta_xip_terms[lfv,j,:,3] = T_ratio**2 *(2*rho[lfv,3,:])
			delta_xip_total[lfv,j,:] = delta_xip_terms[lfv,j,:,0]+delta_xip_terms[lfv,j,:,1]+delta_xip_terms[lfv,j,:,2]+delta_xip_terms[lfv,j,:,3]
			#delta_xip[lfv,j,:] = 2*xip_theory[j,:]*T_ratio*deltaT_ratio + T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+2*rho[lfv,3,:]) 

			err_term1 = T_ratio**2*(rho_err[lfv,0,:]**2+rho_err[lfv,2,:]**2+rho_err[lfv,3,:]**2)
			err_term2 = 0. #T_ratio*alphas[lfv,j]*(rho_err[lfv,1,:]**2+rho_err[lfv,4,:]**2)
			err_delta_xip[lfv,j,:] = ( err_term1 - err_term2 )**0.5
	return delta_xip_total, err_delta_xip, delta_xip_terms

def Calc_delta_xip_cterms():
	# This reads in the <cc> correlation function, where c is the additive shear correction
	# This was produced using the codes, PSFRES_CORRMAP/create_c12_mock.py and Calc_2pt_Stats/calc_xi_w_treecorr.py
	# on TreeCorr.
	delta_xip = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
	err_delta_xip = np.zeros_like( delta_xip )
	for lfv in range(len(LFver)):
		for j in range(num_zbins_tot):
			delta_xip[lfv,j,:], err_delta_xip[lfv,j,:]  = np.loadtxt('LFver%s/delta_xi_sys_LFver%s.dat' %(LFver[lfv],LFver[lfv]), usecols=(3,7), unpack=True) 
	return delta_xip, err_delta_xip



def Set_Mandelbaum_Constraints(alpha):
	NLOS_Cov = 1250
	cosmol_Cov = 'fid'

	scale = np.where(theta_data < 72)[0]
	SNR = np.dot( np.transpose(xip_data[scale]), 
			np.dot( np.linalg.inv(Cov_Mat_uc_Survey[scale[0]:scale[-1]+1,scale[0]:scale[-1]+1]), xip_data[scale] ))		

	# ! EITHER USE A THEORY VECTOR FOR REQUIREMENTS, OR USE SMOOTHED DATA VECTOR !
	if Which_Data == 'Data':
		# Lowest Z_B bin has a wee bit of noise in it. 
		# Smoothing it a wee bit for plotting purposes.
		Sm_xip = gaussian_filter( xip_data, Smooth_Scale )
	elif Which_Data == 'Theory':
		theta_theory, xip_theory_tmp = np.loadtxt('%s/xi_p_smith03revised_zKV450_ZBcut%s' %(data_dir, data_ZBlabel), usecols=(0,1), unpack=True)
		xip_theory = np.interp( theta_data, theta_theory, xip_theory_tmp )
		Sm_xip = np.copy( xip_theory )

	# IF YOU USE 0.1-0.3 DATA/COV THE REQ IS LARGER ON SMALL SCALES AND SMALLER ON LARGE SCALES RELATIVE TO IF WE USED (0.7-0.9)
	# HERE'S THE RATIO: array([ 3.74749516,  2.30591263,  3.22123142,  2.4665902 ,  0.55896293, -0.41969508,  0.39032223,  0.07134366,  0.79189828])
	# SO IN SHORT. NOT OBVIOUS WHICH ONE ZB TO USE. SO GO WITH EITHER. 
	Requirement_134 = T_ratio**(-2) * Sm_xip / (2*SNR)
	Requirement_25  = T_ratio**(-1) * Sm_xip / (2*SNR*alpha)
	return Requirement_134, Requirement_25


def Set_Heymans_Constraints(rho):
	Cov_Mat_uc_Survey_All = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[0:135, 0:135] * Linc_Rescale
	lfv = 0
	delta_xip, _ = Calc_delta_xip_J16( alpha, T_ratio, rho )
	delta_xip = delta_xip[lfv]

	delta_chi2 = np.dot( np.transpose(delta_xip), 
			np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), delta_xip ))		
	return delta_chi2


# Set the requirements on rho_1,2,4 and rho_2,5
if Requirement == "M18":
	Req_134, Req_25 = Set_Mandelbaum_Constraints(alpha.max())  # Just set single alpha val to 0.03 for now.
                                						# need to edit Mandelbaum_* func to work with alpha per bin
elif Requirement == "Z18":
	Req_134 = gaussian_filter( xip_data, Smooth_Scale ) / 10.
	Req_25  = gaussian_filter( xip_data, Smooth_Scale ) / 10.
else:
	print("Requirement must be set to M18 or Z18")
	sys.exit()


def Set_Scales(ax):
	ax.set_xscale('log')
	ax.set_xlim([0.5,300.])
	#ax.set_yscale('log')
	return



# This function plots all 5 rho stats & \delta\xi_+ in separate panels
# with symlog y-axes. A yellow band illustrates the Mandelbaum/Zuntz requirements
# (whichever one is specified at the top)
def Plot_5_Symlog(rho, rho_err, pm):
	
	#delta_xip, err_delta_xip = Calc_delta_xip_J16( alpha, T_ratio, rho, rho_err )
	#delta_xip, err_delta_xip, _ = Calc_delta_xip_H20( T_ratio, rho, rho_err )
	delta_xip, err_delta_xip = Calc_delta_xip_cterms( )

	Deltaxip = True
	if Deltaxip:
		scrollthrough=6
	else:
		scrollthrough=5

	fig = plt.figure(figsize = (24,10))
	gs1 = gridspec.GridSpec(3, 2)
	colors = [ 'magenta', 'darkblue', 'dimgrey', 'orange', 'lawngreen', 'cyan' ] 	
	for i in range(scrollthrough):
			ax = plt.subplot(gs1[i], adjustable='box')
			Set_Scales(ax)

			# Set the requirement bands - rho_2,5 have different requirements than rho_1,3,4
			if i==1 or i==4:
				# rho 2 and 5 have one requirement
				symlogscale = 1e-6
				ax.fill_between(theta_data[:], y1=abs(Req_25)*-1, y2=abs(Req_25)*1, facecolor='yellow') 
				ax.set_ylim([ abs(Req_25[0])*-1, abs(Req_25[0])*1 ])
			elif i==0 or i==2 or i==3:
				symlogscale = 1e-8
				# rho 1,3 and 4 have another requirement
				ax.fill_between(theta_data[:], y1=abs(Req_134)*-1, y2=abs(Req_134)*1, facecolor='yellow') 
				ax.set_ylim([ abs(Req_134[0])*-1, abs(Req_134[0])*1 ])
			else:
				symlogscale = 1e-7
				Req = np.sqrt( np.diagonal(Cov_Mat_uc_Survey) ) / 2.
				ax.fill_between(theta_data[:], y1=abs(Req)*-1, y2=abs(Req)*1, facecolor='yellow') 
				
				
			ax.set_yscale('symlog', linthreshy=symlogscale )
			ax.plot( [0.5,300.], [symlogscale, symlogscale],  'k--' )
			ax.plot( [0.5,300.], [-1*symlogscale, -1*symlogscale], 'k--' )

			# Set the theta bands
			if i ==5:
				ax.set_ylabel(r'$\delta\xi_%s(\theta)$'%pm)
				ax.set_xlabel(r'$\theta$ [arcmin]')
				for lfv in range(len(LFver)):
					idx_max_alpha = np.argmax(alpha[lfv])
					ax.errorbar(theta, delta_xip[lfv,idx_max_alpha,:], yerr=err_delta_xip[lfv,idx_max_alpha,:], 
													color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])
			else:
				ax.set_ylabel(r'$\rho_{%s}(\theta)$'%(i+1))
				if i!=3 and i!=4:
					ax.set_xticks([])
				else:
					ax.set_xlabel(r'$\theta$ [arcmin]')
					# 1, 3, 4 have another requirement
				for lfv in range(len(LFver)):
					ax.errorbar(theta, rho[lfv,i,:], yerr=rho_err[lfv,i,:], color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])

	#fig.suptitle('Lensfit v'+LFver[0])
	if Deltaxip == False:
		plt.subplot(gs1[-1], adjustable='box').set_visible(False)
	plt.subplots_adjust(hspace=0)
	plt.savefig('LFver%s/rho1/Plot_rho%s_CovPatches%sx%s_Require%s_%sxip_Symlog.png'%(LFver[0],pm,Res,Res,Requirement,Which_Data))
	plt.show()
	return
#Plot_5_Symlog(rhop_mean, rhop_err, '+')

# This plots all 5 rho stats on one panel, \delta\xi_+ on a separate panel
# with a broad band illustrating the Mandelbaum/Zuntz requirements
# (whichever one is specified at the top).
def Plot_2Panel(rho, rho_err, pm, lfv):

	#delta_xip, err_delta_xip = Calc_delta_xip_J16( alpha, T_ratio, rho, rho_err )
	delta_xip, err_delta_xip = Calc_delta_xip_cterms( )

	# Just extract the correct LensFit version
	delta_xip = delta_xip[lfv]
	err_delta_xip = err_delta_xip[lfv]

	fig = plt.figure(figsize = (10,10))
	gs1 = gridspec.GridSpec(2, 1)
	colors = [ 'magenta', 'darkblue', 'red', 'orange', 'lawngreen', 'cyan' ] 	
	for i in range(2):
		ax = plt.subplot(gs1[i], adjustable='box')
		Set_Scales(ax)
		ax.set_yscale('log')

		ax.fill_between(theta_data[:], y1=-1*abs(Req_25), y2=abs(Req_25), facecolor='pink') 
		ax.fill_between(theta_data[:], y1=-1*abs(Req_134), y2=abs(Req_134), facecolor='yellow') 
		if i == 0:
			ax.set_ylim( [1.1e-10,1e-4] )
			ax.set_ylabel(r'$|\rho(\theta)|$')
			ax.set_xticks([])
			for j in range(5):
				ax.errorbar(theta, abs(rho[lfv,j,:]), yerr=abs(rho_err[lfv,j,:]), color=colors[j], linewidth=3, label=r'$\rho_%s$'%(j+1) )	
		ax.legend(loc='upper right', frameon=False, ncol=3)
		if i == 1:
			ax.set_ylim( [1.1e-10,1e-4] )
			ax.set_ylabel(r'$\delta\xi_%s(\theta)$'%pm)
			for j in range(num_zbins_tot):
				ax.errorbar(theta, abs(delta_xip[j,:]), yerr=abs(err_delta_xip[j,:]), linewidth=3)
			ax.set_xlabel(r'$\theta$ [arcmin]')

	plt.subplots_adjust(hspace=0, wspace=0)
	plt.savefig('LFver%s/rho1/Plot_Overall-rho%s_CovPatches%sx%s_Require%s_%sxip.png'%(LFver[lfv],pm,Res,Res, Requirement, Which_Data))
	#plt.show()
	return

#Plot_2Panel(rhop_mean, rhop_err, '+', -1)


# This plots the various ingredients of 
def Plot_deltaxips_Only():
	lfv = 0
	zbin= -1 # Plot the delta_xip for this z-bin alone.	

	delta_xip_H, err_delta_xip_H, delta_xip_terms_H = Calc_delta_xip_H20( T_ratio, rhop_mean, rhop_err )
	delta_xip_c, err_delta_xip_c = Calc_delta_xip_cterms( )


	fig = plt.figure(figsize = (10,6))
	gs1 = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs1[0], adjustable='box')
	Set_Scales(ax)
	ax.set_ylabel(r'Components of $\delta\xi_+^{\rm sys}$')
	ax.set_xlabel(r'$\theta$ [arcmin]')

	symlogscale=1e-8
	ax.set_yscale('symlog', linthreshy=symlogscale )
	ax.plot( [0.5,300.], [symlogscale, symlogscale],  'k-' )
	ax.plot( [0.5,300.], [-1*symlogscale, -1*symlogscale], 'k-' )
	ax.set_ylim([-1*1e-5, 1e-5])


	# Plot the diagonal covariance
	Req = np.sqrt( np.diagonal(Cov_Mat_uc_Survey) ) / 2.
	ax.fill_between(theta_data[:], y1=abs(Req)*-1, y2=abs(Req)*1, facecolor='yellow') 

	# Plot the individual ingreidents of the delta_xip in Catherine's model
	ax.plot( theta, delta_xip_terms_H[lfv,zbin,:,0], color='blue', linewidth=2, linestyle=':', 
			label=r'$ 2 \left[{\frac{\overline{\delta T_{\rm PSF}}}{T_{\rm gal}}}\right] \left< e_{\rm obs}^{\rm perfect} e_{\rm obs}^{\rm perfect} \right>$' )
	ax.plot( theta, delta_xip_terms_H[lfv,zbin,:,1], color='blue', linewidth=2, linestyle='--',
			label=r'$\overline{T_{\rm gal}^{-2}} \left< (\delta e_{\rm PSF} \, T_{\rm PSF}) \,  (\delta e_{\rm PSF} \,T_{\rm PSF}) \right>$' )
	ax.plot( theta, delta_xip_terms_H[lfv,zbin,:,2], color='blue', linewidth=2, linestyle='-.',
			label=r'$\overline{T_{\rm gal}^{-2}}\left< (e_{\rm PSF} \, \delta T_{\rm PSF}) \,  (e_{\rm PSF} \, \delta T_{\rm PSF}) \right>$' )
	ax.plot( theta, delta_xip_terms_H[lfv,zbin,:,3], color='blue', linewidth=2, linestyle='-',
			label=r'$2 \overline{T_{\rm gal}^{-2}} \left< (e_{\rm PSF} \, \delta T_{\rm PSF}) \,  (\delta e_{\rm PSF} \, T_{\rm PSF}) \right>$' )

	ax.plot( theta, delta_xip_H[lfv,zbin,:], color='red', linewidth=2, label=r'$\delta\xi_+^{\rm sys}$' ) 
	ax.plot( theta, delta_xip_c[lfv,zbin,:], color='dimgrey', linewidth=2, label=r'$\langle cc \rangle$' )

	ax.legend(loc='upper right', frameon=False, ncol=2)
	plt.subplots_adjust(hspace=0)
	#plt.savefig('LFver%s/rho1/Plot_deltarho%s_CovPatches%sx%s.png'%(LFver[0],pm,Res,Res))
	plt.show()
	return
Plot_deltaxips_Only()


# This plots 1 panel showing the fractional size of the \delta\xi_+
# relative to a noisy theoretical prediction from Marika
def Plot_xip_Plus_rho(rho, rho_err, pm):
	
	fig = plt.figure(figsize = (10,7.5))
	gs1 = gridspec.GridSpec(1, 1)
	# Read in a noisy mock data vector from Marika
	xip_Mock = np.loadtxt('Marika_Noisy_Data_Mock/xipm_nBins1_0.50-300.00_real0')[:7] # Only first 7 points are for xi+

	colors = [ 'magenta', 'darkblue', 'dimgrey', 'orange', 'lawngreen', 'cyan' ] 	

	ax = plt.subplot(gs1[0], adjustable='box')
	Set_Scales(ax)
	ax.set_ylabel(r'$(\xi_+ + \delta_{\xi_+})/\xi_+)$')

	ax.set_xlabel(r'$\theta$ [arcmin]')
	
	for lfv in range(len(LFver)):
		# This delta_xip is eqn 3.9 in Zuntz.
		delta_xip = T_ratio**2*(rho[lfv,0,:7]+rho[lfv,2,:7]+rho[lfv,3,:7]) - T_ratio*alpha*(rho[lfv,1,:7]+rho[lfv,4,:7])
		ax.errorbar(theta[:7], (xip_Mock+delta_xip)/xip_Mock, yerr=rho_err[lfv,0,:7]/xip_Mock, color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])
	ax.plot([0.001,1000], [1,1], 'k--')

	ax.legend(loc='lower right', frameon=False)
	plt.subplots_adjust(hspace=0)
	plt.title(r'%s PSF leakage'%alpha)
	plt.savefig('LFver%s/rho1/Plot_deltarho%s_CovPatches%sx%s.png'%(LFver[0],pm,Res,Res))
	#plt.show()

	return
#Plot_xip_Plus_rho(rhop_mean, rhop_err, '+')





t1 = time.time()
# This function will calculate chi^2 values for many noise realisations under 4 hypotheses:
# null : the chi^2 values obtained for the correct cosmology given no systematics
# sys: the chi^2 values obtained for the correct cosmology given the systematic bias \delta\xi_+
# hi/lo: the chi^2 values obtained for a cosmology that is ~0.004 higher in S_8 than the truth.
# This funcion evaluates if the shift (null-->sys) is subdominant to the shift (null-->hi/lo)
# Turns out it is, therefor systematic is subdominant to this tiny change in cosmology.

def Investigate_chi2(rho):
	Cov_Mat_uc_Survey_All = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[0:135,0:135] * Linc_Rescale
	# [0:135,0:135] pulls out only the xi+ elements.
	lfv = 0
	#delta_xip, _ = Calc_delta_xip_J16( alpha, T_ratio, rho, rho )
	#delta_xip, _ = Calc_delta_xip_H20( T_ratio, rho, rho )
	delta_xip, err_delta_xip = Calc_delta_xip_cterms( )

	delta_xip = np.ndarray.flatten( delta_xip[lfv] )


	# Assemble the theory vector - used to guage signif. of measuring genuine signal.
	# And deviations in this signal from those with high/low values of S_8
	xip_theory_stack_hi = Read_In_Theory_Vector('high')		# High S_8
	xip_theory_stack_lo = Read_In_Theory_Vector('low')		# Low S_8
	xip_theory_stack_fid = Read_In_Theory_Vector('fid')		# Fiducial S_8


	n_noise = 5000
	chi2_null = np.empty( n_noise )		# chi^2 of the null hypothesis (measurement is all noise) for each noise realisation
	chi2_sys = np.empty( n_noise )		# same for the hypothesis that measurement is contaminated by systematic
	chi2_hi = np.empty( n_noise )		# same for the hypothesis that measurement is high S_8
	chi2_lo = np.empty( n_noise )		# same for the hypothesis that measurement is low S_8
	for i in range(n_noise):
		noise = multi_norm.rvs(mean=np.zeros(135), cov=Cov_Mat_uc_Survey_All)
		# chi2 for null hypothesis 
		chi2_null[i] = np.dot( np.transpose(noise), np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), noise ))
		# chi2 for systematic hypothesis
		chi2_sys[i] = np.dot( np.transpose(noise+delta_xip), 
								np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), noise+delta_xip ))
		# chi2 for high S_8 cosmology
		chi2_hi[i] = np.dot( np.transpose(noise+xip_theory_stack_hi-xip_theory_stack_fid), 
								np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), noise+xip_theory_stack_hi-xip_theory_stack_fid ))
		# chi2 for low S_8 cosmology
		chi2_lo[i] = np.dot( np.transpose(noise+xip_theory_stack_lo-xip_theory_stack_fid), 
								np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), noise+xip_theory_stack_lo-xip_theory_stack_fid ))

	# Histogram the chi^2
	def Hist_n_Normalise(chi2, nbins):
		histo_chi2, tmp_bins = np.histogram(chi2, nbins)
		bins_chi2 = tmp_bins[:-1] + (tmp_bins[1] - tmp_bins[0])/2.			# get bin centres
		histo_chi2 = histo_chi2 / simps( histo_chi2, bins_chi2)				# Normalise to a PDF
		mean_chi2 = np.mean( chi2 )
		return histo_chi2, bins_chi2, mean_chi2
	histo_chi2_null, bins_chi2_null, mean_chi2_null = Hist_n_Normalise(chi2_null, 100)
	histo_chi2_sys, bins_chi2_sys, mean_chi2_sys = Hist_n_Normalise(chi2_sys, 100)
	histo_chi2_hi, bins_chi2_hi, mean_chi2_hi = Hist_n_Normalise(chi2_hi, 100)
	histo_chi2_lo, bins_chi2_lo, mean_chi2_lo = Hist_n_Normalise(chi2_lo, 100)
	print("Shift in the mean of the chi2 distributions between null and the following hypotheses is...:" )
	print("sys: ", abs(mean_chi2_null-mean_chi2_sys))
	print("high S8: ", abs(mean_chi2_null-mean_chi2_hi))
	print("low S8: ", abs(mean_chi2_null-mean_chi2_lo))
	


	# Plot the chi^2 distribution
	f, ((ax1)) = plt.subplots(1, 1, figsize=(10,9))
	ax1.plot(bins_chi2_hi, histo_chi2_hi, color='magenta', linewidth=3)
	ax1.plot(bins_chi2_sys, histo_chi2_sys, color='orange', linewidth=3)
	ax1.plot(bins_chi2_null, histo_chi2_null, color='cyan', linewidth=3)

	#ax1.plot( chi2_array, pdf_bf, color='blue', linestyle=':', linewidth=3, label=r'$\rm{DOF}=%.1f \pm %.1f$'%(dof_bf,np.sqrt(dof_err)) )

	# Mark the mean of the chi2 distributions to see how distinguishable they are.
	ax1.plot( [mean_chi2_hi,mean_chi2_hi], [0., histo_chi2_null.max()], 
				color='magenta', linestyle='--', linewidth=3, label=r'$\langle \chi^2_{\Delta S_8} \rangle$' )
	ax1.plot( [mean_chi2_sys,mean_chi2_sys], [0., histo_chi2_null.max()], 
				color='orange', linestyle='--', linewidth=3, label=r'$\langle \chi^2_{\rm sys} \rangle$' )
	ax1.plot( [mean_chi2_null,mean_chi2_null], [0., histo_chi2_null.max()], 
				color='cyan', linestyle='--', linewidth=3, label=r'$\langle \chi^2_{\rm null} \rangle$' )

	#ax1.set_xlim([ chi2_null.min(), chi2_null.max() ])
	ax1.set_xlabel(r'$\chi^2$')
	ax1.set_ylabel(r'PDF$(\chi^2)$')
	ax1.legend(loc='best', frameon=True)
	#plt.savefig('LFver%s/rho1/Plot_PDF-chi2.png'%(LFver[0]))
	plt.show()


	# The following bits of code are no longer used - they compare the chi^2 distributions
	# via a P-value analysis. Instead of doing this, we now just look at the means of the chi^2 distributions
	# (done above). 
	# -------------------------------------------------------------------------------------------- #

	# Fit for the effective DoF of the chi^2 distribution
	def chi2pdf_model(chi2_array, *dof):
		chi2_pdf = chi2.pdf( chi2_array, dof ) / simps( chi2.pdf( chi2_array, dof ), chi2_array)
		# Interpolate the chi^2 to bins of histogram.
		return np.interp( bins_chi2_null, chi2_array, chi2_pdf )	

	chi2_array = np.linspace( 0., 400., 1000 )		# Arbitrary array spanning a suitable range of chi2 values
	p0 = [100]										# Arbitrary first guess at the DoF
	dof_bf, dof_err = curve_fit(chi2pdf_model, 
					chi2_array, histo_chi2_null, p0=p0)	# Best fit dof.
	pdf_bf = chi2.pdf( chi2_array, dof_bf ) / simps( chi2.pdf( chi2_array, dof_bf ), chi2_array)


	# Get the P-values for the chi2_null, chi2_sys, chi2_hi/lo hypotheses, for each noise realisation
	P_null = chi2.cdf( chi2_null, dof_bf )
	P_sys = chi2.cdf( chi2_sys, dof_bf )	
	P_hi = chi2.cdf( chi2_hi, dof_bf )	
	P_lo = chi2.cdf( chi2_lo, dof_bf )	


	def Plot_Pvalue_Histo(P, label):
		# Plot the distribution of P-values for the noise-realisations
		if len(P.shape) == 1:
			P = np.reshape(P, (1,len(P)))

		colours = ['magenta', 'dimgrey', 'orange', 'cyan', 'darkblue']		
		f, ((ax1)) = plt.subplots(1, 1, figsize=(10,9))
		dummy, bin_edges = np.histogram(P[0,:], 100)		# Use this just to get bin edges used for all distributions
		bins_P = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2.
		histo_P = np.zeros([ P.shape[0], len(bins_P) ])
		cumul_P = np.zeros([ P.shape[0], len(bins_P) ])
		for i in range(P.shape[0]):
			histo_P[i,:] = np.histogram(P[i,:], bin_edges)[0]
			cumul_P[i,:] = np.cumsum(histo_P[i,:])
			ax1.bar(bins_P, histo_P[i,:], width=(bins_P[1]-bins_P[0]), color=colours[i], 
					edgecolor=colours[i], alpha=0.5, label=label[i])
			#ax1.step(bins_P, cumul_P[i,:], color=colours[i], label=label[i])

		ax1.set_xlabel(r'$P$') # -P_{\rm null}
		ax1.set_ylabel(r'PDF$(P)$')
		#ax1.set_ylabel(r'CDF$(P)$')		
		ax1.legend(loc='best', frameon=False)
		#plt.savefig('LFver%s/rho1/Plot_PDF-Pvalues_%s.png'%(LFver[0], label))
		#plt.show()
		return histo_P, cumul_P
	P_stack = np.vstack(( P_hi, P_lo, P_sys, P_null ))
	histo_P_stack, cumul_P_stack = Plot_Pvalue_Histo(P_stack, [r'high $S_8$', r'low $S_8$', 'sys', 'null'] )

	# Do KS tests between the P-values of the chi^2's of the hi/lo/null and systematic hypotheses
	# Could compare the PDFs or the CDFs. Turns out it makes a big difference
	# with PDFs varying wildly with the noise realisations, and the CDFs 
	# being completely indistinguishable due to the ~flatness of the PDFs.
	# This is part of the reason we no longer do this P-value analysis.
	# "-------------------PDF KS Values---------------------" 
	KS_sys_null = ks_2samp(histo_P_stack[2,:], histo_P_stack[3,:])
	KS_hi_null = ks_2samp(histo_P_stack[0,:], histo_P_stack[3,:])
	KS_lo_null = ks_2samp(histo_P_stack[1,:], histo_P_stack[3,:])
	# "-------------------CDF KS Values---------------------" 
	#KS_sys_null = ks_2samp(cumul_P_stack[2,:], cumul_P_stack[3,:])
	#KS_hi_null = ks_2samp(cumul_P_stack[0,:], cumul_P_stack[3,:])
	#KS_lo_null = ks_2samp(cumul_P_stack[1,:], cumul_P_stack[3,:])
	# "------------------------------------------------------"

	# -------------------------------------------------------------------------------------------- #


	return histo_P_stack, KS_sys_null, KS_hi_null, KS_lo_null
histo_P_stack, KS_sys_null, KS_hi_null, KS_lo_null = Investigate_chi2(rhop_mean)
t2 = time.time()
print(" It took %.0f seconds." %(t2-t1))



