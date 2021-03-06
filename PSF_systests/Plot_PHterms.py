# 14/05/2020: B. M. Giblin, Edinburgh
# Plot the additive & multiplicative systematics to \xi+ as predicted in the
# Pauli-Henriksson model detailed in eqn 10 of Giblin, Heymans, Asgari et al. (2020).

#import sys
#sys.path.insert(0, './')
from Functions_4_Plotting_PSFstats import Calc_Important_Tquantities, Read_rho_Or_PH, Select_Patch, Read_In_Theory_Vector

import numpy as np
import pylab as plt
from astropy.io import fits
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import chi2
from scipy.stats import multivariate_normal as multi_norm

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
        'size'   : 26} # 12

plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

# Which Lensfit version catalogue were the PHterms calculated for?
LFver = ["321"] #["309b",  "319b", "319c", "319d", "321"] # "319",
# boundaries of the KiDS-1000 cosmic shear redshift bins.
zbounds = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2 ]
num_zbins = len(zbounds)-1
num_zbins_tot = np.sum( range(num_zbins+1) )    # Number source bins including cross-bins

# "T-quantities" here refers to the (delta)T_(PSF/gal) estimates - the PSF sizes defined in Giblin et al. (2020).
# Arrays to store T-quantities (2nd dimensions is len 2,
# where 0=mean, 1=error)
deltaT_ratio = np.zeros([ len(LFver), 2, num_zbins ])          # deltaT_PSF/T_gal for the 5 zbins
Tg_invsq = np.zeros_like(deltaT_ratio)                         # 1/T_gal^2 for the 5 zbins
deltaT_ratio_tot = np.zeros([ len(LFver), 2, num_zbins_tot ])  # same as above for all auto & cross bins
Tg_invsq_tot = np.zeros_like(deltaT_ratio_tot)                 # " "

Read_Tvalues = True # If True, read in pre-saved values of gal & PSF size quantities 
if Read_Tvalues:
        print('Reading in pre-saved T-quantities for each ZB bin.')
        for lfv in range(len(LFver)):
                deltaT_ratio_tot[lfv,0,:], deltaT_ratio_tot[lfv,1,:],Tg_invsq_tot[lfv,0,:], Tg_invsq_tot[lfv,1,:] = np.loadtxt('LFver%s/PHterms/deltaTratio_TgInvSq_MeansAndErrors_%stomobins.dat' %(LFver[lfv],num_zbins_tot),usecols=(0,1,2,3), unpack=True)
                
else:
        # Calculate them
        print('Calculating T-quantities for each ZB bin')
        nboot = 30
        for lfv in range(len(LFver)):
                # The following function avg's the T_gal,T_PSF & (T_PSF-T_gal) across the survey and redshifts:
                deltaT_ratio[lfv], Tg_invsq[lfv], deltaT_ratio_tot[lfv], Tg_invsq_tot[lfv] = Calc_Important_Tquantities(LFver[lfv],zbounds, nboot)
                np.savetxt('LFver%s/PHterms/deltaTratio_TgInvSq_MeansAndErrors_%stomobins.dat' %(LFver[lfv],num_zbins_tot),
                           np.c_[ deltaT_ratio_tot[lfv,0,:], deltaT_ratio_tot[lfv,1,:],
                                  Tg_invsq_tot[lfv,0,:], Tg_invsq_tot[lfv,1,:] ],
                           header='# < deltaT_PSF / T_gal >, SIGMA[ deltaT_PSF / T_gal ], < T_gal^-2 > , SIGMA[ T_gal^-2  ]')

# alpha is the PSF leakage defined in eqn 3.17 of Jarvis et al. (2016).
# In this code, is only used in the Jarvis prescription for the systematic to xi+.
# Read in the alpha values for each shear component and tomo-bin:
Use_alpha_per_bin = True
if Use_alpha_per_bin:
        from Functions_4_Plotting_PSFstats import Read_alpha_per_bin
        alpha, alpha_err = Read_alpha_per_bin(["321", "309c"])   # Read in alpha values for both LFver
else:
        alpha = np.zeros([ len(LFver), num_zbins_tot ])+0.03     # worse case scenario PSF leakage is ~0.03
        alpha_err=np.zeros([ len(LFver), num_zbins_tot ])

ThBins = 9  # number of theta bins used by TreeCorr in PHterms.
Res = 7     # Res*Res is the number of patches used in jackknife realisations.
# Read in the mean and error on the PH stats + the labels for plotting.
Plot_Labels, theta, php_mean, php_err = Read_rho_Or_PH(LFver, 'ph', ThBins, Res) 
#sys.exit()


def Marika_Cov(mCov):
        # The KiDS-1000 cosmic shear covariance used in Asgari et al. (2020)
        # if mCov is True, reads covariance include m-uncertainty
        # else reads cov excluding this.
        if mCov:
                keyword = 'with'
        else:
                keyword = 'no'
        cov_fname = 'Marikas_CovMat/xipm_KIDS1000_BlindA_%s_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits' %keyword
        f = fits.open(cov_fname)
        theta_cov = f[2].data['ANG'][0:ThBins] # The theta values the cov is defined at (in arcmin)
        cov = f[1].data
        return theta_cov, cov #[81:90, 81:90]




# THE FOLLOWING FUNCTIONS ARE DIFFERENT PRESCRIPTIONS FOR THE SYSTEMATIC TO THE XI+ SHEAR CORRELATION FUNCTION: delta_xip


def Calc_delta_xip_J16(alphas, T_ratio, rho, rho_err): 		# !!! Jarvis (2016) expression for PSF systematic
	xip_theory = Read_In_Theory_Vector('fid', theta_data)	# fiducial cosmology xi+ theory vector
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



def Calc_delta_xip_H20(ph, ph_err): 		# !!! The Paulin-Henriksson model appearing in Giblin (2020) eqn 10.
                                                # Credit to Catherine Heymans for its derivation
        theta_theory, xip_theory = Read_In_Theory_Vector('fid')	
        xip_theory = np.reshape(xip_theory, (num_zbins_tot,len(theta_theory))) # Reshape to be [num_zbins,ntheta]

        # Calculate the additive shear bias for each lensfit version and in each tomographic bin 
        delta_xip_total = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
        delta_xip_terms = np.zeros([ len(LFver), num_zbins_tot, ThBins, 4 ]) # Store the separate ingredients of the total delta_xip

        err_delta_xip = np.zeros_like( delta_xip_total )
        err_delta_xip_terms = np.zeros_like( delta_xip_terms )
        for lfv in range(len(LFver)):
                for j in range(num_zbins_tot):
                        # The 4 individual terms in \delta\xi+ (eqn 10. Giblin et al. 2020)
                        # Remember, in deltaT_ratio_tot[lfv,i,:], i=0 for mean, i=1 for error
                        # same goes for Tg_invsq_tot.
                        # Also interpolate the theory vector onto the theta bins of the PH-stats
                        # (should only be a minor change in the theta values).
                        delta_xip_terms[lfv,j,:,0] = 2* np.interp(theta, theta_theory, xip_theory[j,:]) *deltaT_ratio_tot[lfv,0,j]
                        delta_xip_terms[lfv,j,:,1] =   Tg_invsq_tot[lfv,0,j] *(ph[lfv,0,:])
                        delta_xip_terms[lfv,j,:,2] = 2*Tg_invsq_tot[lfv,0,j] *(ph[lfv,1,:])
                        delta_xip_terms[lfv,j,:,3] =   Tg_invsq_tot[lfv,0,j] *(ph[lfv,2,:])
                        # The total
                        delta_xip_total[lfv,j,:] = delta_xip_terms[lfv,j,:,0]+delta_xip_terms[lfv,j,:,1]+delta_xip_terms[lfv,j,:,2]+delta_xip_terms[lfv,j,:,3]
                        # Following terms come from error propagation of eqn 10. (Giblin, Heymans et al. 2020)
                        err_delta_xip_terms[lfv,j,:,0] = 2 * np.interp(theta, theta_theory, xip_theory[j,:]) * deltaT_ratio_tot[lfv,1,j]
                        scale = [1,2,1] # PH term 2 has a factor 2 in it, others factor 1.
                        for t in range(1,4): # cycle through 3 ph terms - same error form
                                part1 = scale[t-1] * Tg_invsq_tot[lfv,1,j]**2 * ph[lfv,t-1,:]**2  
                                part2 = scale[t-1] * Tg_invsq_tot[lfv,0,j]**2 * ph_err[lfv,t-1,:]**2
                                err_delta_xip_terms[lfv,j,:,t] = (part1 + part2)**0.5
                        			
                        err_delta_xip[lfv,j,:] = ( err_delta_xip_terms[lfv,j,:,0]**2 + err_delta_xip_terms[lfv,j,:,1]**2 + err_delta_xip_terms[lfv,j,:,2]**2 + err_delta_xip_terms[lfv,j,:,3]**2 )**0.5
        return delta_xip_total, err_delta_xip, delta_xip_terms, err_delta_xip_terms


def Calc_delta_xip_cterms():
	# This reads in the <ePSF(x,y) ePSF(x,y)> correlation function.
	# This was produced using the codes, PSFRES_CORRMAP/create_c12_mock.py and Calc_2pt_Stats/calc_xi_w_treecorr.py
	# on TreeCorr.
	delta_xip = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
	err_delta_xip = np.zeros_like( delta_xip )
	for lfv in range(len(LFver)):
		for j in range(num_zbins_tot):
			theta_dxip, delta_xip[lfv,j,:], err_delta_xip[lfv,j,:]  = np.loadtxt('LFver%s/delta_xi_sys_LFver%s.dat' %(LFver[lfv],LFver[lfv]), usecols=(1,3,7), unpack=True) 
	return theta_dxip, delta_xip, err_delta_xip

def Calc_delta_xip_Bacon():
        # The Bacon et al. (2010) PSF systematic: alpha^2 * xi_+^{PSF,PSF}
        # for the 15 tomographic bins (eqn 13. of Giblin et al. 2020).
        # Note this exists only for LFver svn_309c.
        Bacon_DIR = '/disk09/KIDS/K1000_TWO_PT_STATS/OUTSTATS/CSys/'
        delta_xip = np.zeros([ num_zbins_tot, ThBins ])
        delta_xip_err = np.zeros_like( delta_xip )
        k=0
        for i in range( num_zbins ):
                for j in range( num_zbins ):
                        if i<=j:
                                theta_dxi, delta_xip[k,:], delta_xip_err[k,:] = np.loadtxt('%s/CSys_5Z_%s_%s_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.dat' %(Bacon_DIR,i+1,j+1), usecols=(1,3,5), unpack=True)
                                #delta_xip[k,:] = np.interp( theta, tmp_theta, tmp_dxi )
                                # Don't need to interp onto the theta bins of the PH-stats
                                # Since the theta scales differ by ~<1%.
                                k+=1
        return theta_dxi, delta_xip, delta_xip_err
                                           

# Function that is called just to set the theta limits and scale of plots
def Set_Scales(ax):
	ax.set_xscale('log')
	ax.set_xlim([0.5,350.])
	#ax.set_yscale('log')
	return



# This plots the various ingredients of delta_xip
# (produces Figure 3 of Giblin et al. 2020).
# Input: which redshift bin to produce the plot for?
def Plot_deltaxips_Only(zbin):
        theta_cov, cov = Marika_Cov(True)
        cov = cov[81:90, 81:90]*0.16 # Pulls out final zbin covariance - contact C. Heymans for info on factor 0.16.
        
        lfv = 0
        #zbin= -1 # Plot the delta_xip for this z-bin alone.	

        # Return the various delta_xip prescriptions
        delta_xip_H, err_delta_xip_H, delta_xip_terms_H, err_delta_xip_terms_H = Calc_delta_xip_H20( php_mean, php_err )
        theta_dxip_c, delta_xip_c, err_delta_xip_c = Calc_delta_xip_cterms( )
        theta_dxip_B, delta_xip_B, err_delta_xip_B = Calc_delta_xip_Bacon()

        fig = plt.figure(figsize = (11,10)) #(8,6)
        gs1 = gridspec.GridSpec(2, 1, height_ratios=[0.05,0.95]) #[0.15,0.85]
        
        ax = plt.subplot(gs1[1], adjustable='box')
        Set_Scales(ax)
        ax.set_ylabel(r'Components of $\delta\xi_+^{\rm sys}$')
        ax.set_xlabel(r'$\theta$ [arcmin]')

        symlogscale=1e-10
        ax.set_yscale('symlog', linthreshy=symlogscale )
        ax.plot( [0.5,400.], [symlogscale, symlogscale],  'k-' )
        ax.plot( [0.5,400.], [-1*symlogscale, -1*symlogscale], 'k-' )
        ax.set_ylim([-1*1e-5, 1e-5])


        # Plot the diagonal covariance
        Req = np.sqrt( np.diagonal(cov) ) / 2.
        ax.fill_between(theta_cov, y1=abs(Req)*-1, y2=abs(Req)*1, facecolor='yellow') 

	# Plot the individual ingreidents of the delta_xip in the PH model
        ax.errorbar( 10**(np.log10(theta)-0.05), delta_xip_terms_H[lfv,zbin,:,0],
                     yerr=err_delta_xip_terms_H[lfv,zbin,:,0],
                     color='lime', linewidth=2, linestyle=':', 
		     label=r'$ 2 \left[{\frac{\overline{\delta T_{\rm PSF}}}{T_{\rm gal}}}\right] \left< e_{\rm obs}^{\rm perfect} e_{\rm obs}^{\rm perfect} \right>$' )
        
        ax.errorbar( theta, delta_xip_terms_H[lfv,zbin,:,1],
                     yerr=err_delta_xip_terms_H[lfv,zbin,:,1],
                     color='green', linewidth=2, linestyle='-.',
		     label=r'$\left[ \,\overline{\frac{1}{T_{\rm gal}}}\,\right]^2 \left< (e_{\rm PSF} \, \delta T_{\rm PSF}) \,  (e_{\rm PSF} \, \delta T_{\rm PSF}) \right>$' )
        
        ax.errorbar( 10**(np.log10(theta)+0.025), delta_xip_terms_H[lfv,zbin,:,2],
                     yerr=err_delta_xip_terms_H[lfv,zbin,:,2],
                     color='cyan', linewidth=2, linestyle='-',
		     label=r'$2 \left[ \,\overline{\frac{1}{T_{\rm gal}}}\,\right]^2  \left< (e_{\rm PSF} \, \delta T_{\rm PSF}) \,  (\delta e_{\rm PSF} \, T_{\rm PSF}) \right>$' )

        ax.errorbar( 10**(np.log10(theta)-0.025), delta_xip_terms_H[lfv,zbin,:,3],
                     yerr=err_delta_xip_terms_H[lfv,zbin,:,3],
                     color='blue', linewidth=2, linestyle='--',
                     label=r'$\left[ \,\overline{\frac{1}{T_{\rm gal}}}\,\right]^2 \left< (\delta e_{\rm PSF} \, T_{\rm PSF}) \,  (\delta e_{\rm PSF} \,T_{\rm PSF}) \right>$' )

        # Plot the overall PH delta_xip systematic (eqn 10 of Giblin+20)
        ax.errorbar( 10**(np.log10(theta)+0.05), delta_xip_H[lfv,zbin,:],
                     yerr=err_delta_xip_H[lfv,zbin,:], color='lawngreen', linewidth=3, #'red', 2
                     label=r'$\delta\xi_+^{\rm sys}$' ) # label=r'PH08')

        # Plot the Bacon et al. systematic (eqn 13 of Giblin+20)
        #ax.errorbar( theta_dxip_B, delta_xip_B[zbin,:], yerr=err_delta_xip_B[zbin,:],
        #             color='orange', linewidth=3, #2
        #             label=r'$\alpha^2 \xi^{\rm PSF,PSF}$' ) #label=r'B03')

        # Plot the flux-dependent PSF residual systematic: (Sec. 3.4.1 of Giblin+20). 
        ax.errorbar( theta_dxip_c, delta_xip_c[lfv,zbin,:], yerr=err_delta_xip_c[lfv,zbin,:],
                     color='magenta', linewidth=3, #2
                     label=r'$\langle \delta\epsilon^{\rm PSF}(x,y) \, \delta\epsilon^{\rm PSF}(x,y) \rangle$' )

        
        ax.legend(bbox_to_anchor=(0.5, 1.18), loc='upper center', frameon=False, ncol=2, fontsize=26)
        # bbox_to_anchor=(0.5, 1.37)
        #fontsize=19)
        # bbox_to_anchor=(0.5, 1.5)
        plt.subplots_adjust(hspace=0)
        #plt.savefig('/home/bengib/Poster_Plot.png')
        plt.savefig('LFver%s/PHterms/Plot_deltaxip_CovPatches%sx%s_zbin%s.png'%(LFver[0],Res,Res,zbin))
        plt.show()
        return
#for i in range(num_zbins_tot):
#        Plot_deltaxips_Only(i)
Plot_deltaxips_Only(num_zbins_tot-1)



# This function will calculate chi^2 values for many noise realisations under 4 hypotheses:
# null : the chi^2 values obtained for the correct cosmology given no systematics
# sys: the chi^2 values obtained for the correct cosmology given the systematic bias \delta\xi_+
# hi/lo: the chi^2 values obtained for a cosmology that is ~sigma_shift*sigma higher/lower in S_8 than the truth.
# where sigma_shift is an input argument; we've tried the following values for sigma_shift:  (0.05,0.1,0.2,0.4).
# This funcion evaluates if the shift (null-->sys) is subdominant to the shift (null-->hi/lo)
# first argument 'rho' is the delta_xip systematic array.

def Investigate_chi2(rho, sigma_shift):
        theta_cov, cov_All = Marika_Cov(False)
        cov_All = cov_All[0:135,0:135]
        # [0:135,0:135] pulls out only the xi+ elements.   
        
        # And deviations in this signal from those with high/low values of S_8
        # 'higher/lower' are a 0.004 change in S_8 (0.2 sigma_k1000)
        # 'high/low' are a 0.002 change in S_8 (0.1 sigma_k1000),
        # 'midhigh/midlow' are a 0.001 change in S_8 (0.05 sigma_k1000) 
        theta_hi, xip_theory_stack_hi = Read_In_Theory_Vector('high%ssigma' %sigma_shift)		# High S_8
        theta_low, xip_theory_stack_lo = Read_In_Theory_Vector('low%ssigma' %sigma_shift)		# Low S_8
        theta_fid, xip_theory_stack_fid = Read_In_Theory_Vector('fid')		# Fiducial S_8

        lfv = 0
        # !!! CHANGE THIS LINES TO PICK WHICH delta_xip PRESCRIPTION YOU WISH TO TEST.
        #delta_xip, _ = Calc_delta_xip_J16( alpha, T_ratio, rho, rho )
        delta_xip,_,_,_ = Calc_delta_xip_H20( rho, rho )
        #theta_dxip, delta_xip,_ = Calc_delta_xip_cterms( )
        delta_xip = delta_xip[lfv]
        #theta_dxip, delta_xip,_ = Calc_delta_xip_Bacon()
        for i in range(num_zbins_tot):
                delta_xip[i,:] = np.interp( theta_fid, theta, delta_xip[i,:] )
        delta_xip = np.ndarray.flatten( delta_xip )
        
        n_noise = 5000
        chi2_null = np.empty( n_noise )		# chi^2 of the null hypothesis (measurement is all noise) for each noise realisation
        chi2_sys = np.empty( n_noise )		# same for the hypothesis that measurement is contaminated by systematic
        chi2_hi = np.empty( n_noise )		# same for the hypothesis that measurement is high S_8
        chi2_lo = np.empty( n_noise )		# same for the hypothesis that measurement is low S_8
        for i in range(n_noise):
                np.random.seed(i)
                noise = multi_norm.rvs(mean=np.zeros(135), cov=cov_All)
                # chi2 for null hypothesis 
                chi2_null[i] = np.dot( np.transpose(noise), np.dot( np.linalg.inv(cov_All), noise ))
                # chi2 for systematic hypothesis
                chi2_sys[i] = np.dot( np.transpose(noise+delta_xip), 
				      np.dot( np.linalg.inv(cov_All), noise+delta_xip ))
                # chi2 for high S_8 cosmology
                chi2_hi[i] = np.dot( np.transpose(noise+xip_theory_stack_hi-xip_theory_stack_fid), 
				     np.dot( np.linalg.inv(cov_All), noise+xip_theory_stack_hi-xip_theory_stack_fid ))
                # chi2 for low S_8 cosmology
                chi2_lo[i] = np.dot( np.transpose(noise+xip_theory_stack_lo-xip_theory_stack_fid), 
				     np.dot( np.linalg.inv(cov_All), noise+xip_theory_stack_lo-xip_theory_stack_fid ))

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
        return abs(mean_chi2_null-mean_chi2_sys), abs(mean_chi2_null-mean_chi2_hi), abs(mean_chi2_null-mean_chi2_lo)
t1 = time.time()
sigma_shifts = [0.10] #, 0.15, 0.20, 0.3, 0.4]
ntrials = 1  # NOTE: Only set this to >1 if np.random.seed is commented out in the function
             # otherwise it will just return the exact same delta chi^2 shifts for every trial.
delta_chi2_sys = np.zeros([len(sigma_shifts), ntrials])
delta_chi2_hi = np.zeros_like(delta_chi2_sys)
delta_chi2_lo =	np.zeros_like(delta_chi2_sys)
# Uncomment these lines to run through multiple sigma_shift values
#for j in range(len(sigma_shifts)):
#        for i in range(ntrials):
#                delta_chi2_sys[j,i], delta_chi2_hi[j,i], delta_chi2_lo[j,i] = Investigate_chi2(php_mean, sigma_shifts[j])
t2 = time.time()
print(" It took %.0f seconds." %(t2-t1))

