# 14/05/2020: B. M. Giblin, Edinburgh
# Plot the additive systematic to \xi+ as predicted in the
# Pauli-Henriksson model detailed in 3.4 of Giblin, Heymans et al. (2020).

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
        'size'   : 14} # 19

plt.rc('font', **font)
plt.rcParams["mathtext.fontset"] = "cm"

LFver = ["321"] #["309b",  "319b", "319c", "319d", "321"] # "319",

zbounds = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2 ]
num_zbins = len(zbounds)-1
num_zbins_tot = np.sum( range(num_zbins+1) )    # Number source bins including cross-bins

# Arrays to store T-quantities (2nd dimensions is len 2,
# where 0=mean, 1=error
deltaT_ratio = np.zeros([ len(LFver), 2, num_zbins ])
Tg_invsq = np.zeros_like(deltaT_ratio)
deltaT_ratio_tot = np.zeros([ len(LFver), 2, num_zbins_tot ])
Tg_invsq_tot = np.zeros_like(deltaT_ratio_tot)

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
                deltaT_ratio[lfv], Tg_invsq[lfv], deltaT_ratio_tot[lfv], Tg_invsq_tot[lfv] = Calc_Important_Tquantities(LFver[lfv],zbounds, nboot)
                np.savetxt('LFver%s/PHterms/deltaTratio_TgInvSq_MeansAndErrors_%stomobins.dat' %(LFver[lfv],num_zbins_tot),
                           np.c_[ deltaT_ratio_tot[lfv,0,:], deltaT_ratio_tot[lfv,1,:],
                                  Tg_invsq_tot[lfv,0,:], Tg_invsq_tot[lfv,1,:] ],
                           header='# < deltaT_PSF / T_gal >, SIGMA[ deltaT_PSF / T_gal ], < T_gal^-2 > , SIGMA[ T_gal^-2  ]')
        
# Read in the alpha values for each shear component and tomo-bin
Use_alpha_per_bin = False
alpha = np.zeros([ len(LFver), num_zbins_tot ])
if Use_alpha_per_bin:
        from Functions_4_Plotting_PSFstats import Read_alpha_per_bin
        Read_alpha_per_bin(LFver)
else:
	alpha += 0.03     # worse case scenario PSF leakage is 0.03;


ThBins = 9
Res = 7
# Read in the mean and error on the PH stats + the labels for plotting.
Plot_Labels, theta, php_mean, php_err = Read_rho_Or_PH(LFver, 'ph', ThBins, Res) 
#sys.exit()

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


                                                        
def Calc_delta_xip_J16(alphas, T_ratio, rho, rho_err): 		# !!! Jarvis (2016) expression for PSF systematic
	xip_theory = Read_In_Theory_Vector('fid', theta_data)	
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



def Calc_delta_xip_H20(ph, ph_err): 		# !!! Heymans' (2020) derivation for PSF systematic
        xip_theory = Read_In_Theory_Vector('fid', theta_data)	
        xip_theory = np.reshape(xip_theory, (num_zbins_tot,len(theta))) # Reshape to be [num_zbins,ntheta]

        # Calculate the additive shear bias for each lensfit version and in each tomographic bin (ie- value of alpha)
        delta_xip_total = np.zeros([ len(LFver), num_zbins_tot, ThBins ])
        delta_xip_terms = np.zeros([ len(LFver), num_zbins_tot, ThBins, 4 ]) # Store the separate ingredients of the total delta_xip

        err_delta_xip = np.zeros_like( delta_xip_total )
        err_delta_xip_terms = np.zeros_like( delta_xip_terms )
        for lfv in range(len(LFver)):
                for j in range(num_zbins_tot):
                        # The 4 individual terms in \delta\xi+ (eqn 10. Giblin et al. 2020)
                        # Remember, in deltaT_ratio_tot[lfv,i,:], i=0 for mean, i=1 for error
                        # same goes for Tg_invsq_tot. 
                        delta_xip_terms[lfv,j,:,0] = 2*xip_theory[j,:]*deltaT_ratio_tot[lfv,0,j]
                        delta_xip_terms[lfv,j,:,1] =   Tg_invsq_tot[lfv,0,j] *(ph[lfv,0,:])
                        delta_xip_terms[lfv,j,:,2] = 2*Tg_invsq_tot[lfv,0,j] *(ph[lfv,1,:])
                        delta_xip_terms[lfv,j,:,3] =   Tg_invsq_tot[lfv,0,j] *(ph[lfv,2,:])
                        # The total
                        delta_xip_total[lfv,j,:] = delta_xip_terms[lfv,j,:,0]+delta_xip_terms[lfv,j,:,1]+delta_xip_terms[lfv,j,:,2]+delta_xip_terms[lfv,j,:,3]
                        # Following terms come from error propagation of eqn 10. (Giblin, Heymans et al. 2020)
                        err_delta_xip_terms[lfv,j,:,0] = 2 * xip_theory[j,:] * deltaT_ratio_tot[lfv,1,j]
                        scale = [1,2,1] # PH term 2 has a factor 2 in it, others factor 1.
                        for t in range(1,4): # cycle through 3 ph terms - same error form
                                part1 = scale[t-1] * Tg_invsq_tot[lfv,1,j]**2 * ph[lfv,t-1,:]**2  
                                part2 = scale[t-1] * Tg_invsq_tot[lfv,0,j]**2 * ph_err[lfv,t-1,:]**2
                                err_delta_xip_terms[lfv,j,:,t] = (part1 + part2)**0.5
                        			
                        err_delta_xip[lfv,j,:] = ( err_delta_xip_terms[lfv,j,:,0]**2 + err_delta_xip_terms[lfv,j,:,1]**2 + err_delta_xip_terms[lfv,j,:,2]**2 + err_delta_xip_terms[lfv,j,:,3]**2 )**0.5
        return delta_xip_total, err_delta_xip, delta_xip_terms, err_delta_xip_terms


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


def Set_Heymans_Constraints(rho):
	Cov_Mat_uc_Survey_All = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[0:135, 0:135] * Linc_Rescale
	lfv = 0
	delta_xip, _ = Calc_delta_xip_J16( alpha, T_ratio, rho )
	delta_xip = delta_xip[lfv]

	delta_chi2 = np.dot( np.transpose(delta_xip), 
			np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), delta_xip ))		
	return delta_chi2



def Set_Scales(ax):
	ax.set_xscale('log')
	ax.set_xlim([0.5,400.])
	#ax.set_yscale('log')
	return



# This plots the various ingredients of 
def Plot_deltaxips_Only(zbin):
        lfv = 0
        #zbin= -1 # Plot the delta_xip for this z-bin alone.	
        
        delta_xip_H, err_delta_xip_H, delta_xip_terms_H, err_delta_xip_terms_H = Calc_delta_xip_H20( php_mean, php_err )
        delta_xip_c, err_delta_xip_c = Calc_delta_xip_cterms( )


        fig = plt.figure(figsize = (10,6))
        gs1 = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs1[0], adjustable='box')
        Set_Scales(ax)
        ax.set_ylabel(r'Components of $\delta\xi_+^{\rm sys}$')
        ax.set_xlabel(r'$\theta$ [arcmin]')

        symlogscale=1e-10
        ax.set_yscale('symlog', linthreshy=symlogscale )
        ax.plot( [0.5,400.], [symlogscale, symlogscale],  'k-' )
        ax.plot( [0.5,400.], [-1*symlogscale, -1*symlogscale], 'k-' )
        ax.set_ylim([-1*1e-5, 1e-1])


        # Plot the diagonal covariance
        Req = np.sqrt( np.diagonal(Cov_Mat_uc_Survey) ) / 2.
        ax.fill_between(theta_data, y1=abs(Req)*-1, y2=abs(Req)*1, facecolor='yellow') 

	# Plot the individual ingreidents of the delta_xip in Catherine's model
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

        ax.errorbar( 10**(np.log10(theta)+0.05), delta_xip_H[lfv,zbin,:],
                     yerr=err_delta_xip_H[lfv,zbin,:], color='red',
                     linewidth=2, label=r'$\delta\xi_+^{\rm sys}$' )
        
        ax.errorbar( theta, delta_xip_c[lfv,zbin,:], yerr=err_delta_xip_c[lfv,zbin,:],
                     color='magenta', linewidth=2,
                     label=r'$\langle \delta\epsilon^{\rm PSF}(x,y) \, \delta\epsilon^{\rm PSF}(x,y) \rangle$' )

        ax.legend(loc='upper right', frameon=False, ncol=2)
        plt.subplots_adjust(hspace=0)
        plt.savefig('LFver%s/PHterms/Plot_deltaxip_CovPatches%sx%s_zbin%s.png'%(LFver[0],Res,Res,zbin))
        plt.show()
        return
#for i in range(num_zbins_tot):
#        Plot_deltaxips_Only(i)
Plot_deltaxips_Only(num_zbins_tot-1)


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
	delta_xip,_,_,_ = Calc_delta_xip_H20( rho, rho )
	#delta_xip, err_delta_xip = Calc_delta_xip_cterms( )

	delta_xip = np.ndarray.flatten( delta_xip[lfv] )


	# Assemble the theory vector - used to guage signif. of measuring genuine signal.
	# And deviations in this signal from those with high/low values of S_8
        # 'high/low' are a 0.002 change in S_8 (0.1 sigma_k1000),
        # 'midhigh/midlow' are a 0.001 change in S_8 (0.05 sigma_k1000) 
	xip_theory_stack_hi = Read_In_Theory_Vector('midhigh', theta_data)		# High S_8
	xip_theory_stack_lo = Read_In_Theory_Vector('midlow', theta_data)		# Low S_8
	xip_theory_stack_fid = Read_In_Theory_Vector('fid', theta_data)		# Fiducial S_8


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
#histo_P_stack, KS_sys_null, KS_hi_null, KS_lo_null = Investigate_chi2(php_mean)
t2 = time.time()
print(" It took %.0f seconds." %(t2-t1))



