import numpy as np
import pylab as plt
from matplotlib import rc
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
alpha=0.03 		# worse case scenario PSF leakage is 0.03;
			   	# but you still need to calculate this term from data
T_ratio = 1. 	# Will ultimately get this column from Lance's updated Lensfit cats [I think?]



LFver = ["321"] #["309b",  "319b", "319c", "319d", "321"] # "319",
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
Cov_inDIR = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/Lincs_CovMat/'
Cov_Mat_uc_Survey = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[126:135, 126:135] * Linc_Rescale		# [0:9, 0:9] This extracts xi+ Cov in lowest bin
																											# [126:135, 126:135] for xi+ in highest bin

def Set_Mandelbaum_Constraints():
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
	delta_xip = T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alpha*(rho[lfv,1,:]+rho[lfv,4,:])
	# append to itself 15 times for the 15 bins (delta_xip same for every bin).
	delta_xip_stack = np.tile(delta_xip, 15)	
	delta_chi2 = np.dot( np.transpose(delta_xip_stack), 
			np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), delta_xip_stack ))		
	return delta_chi2


# Set the requirements on rho_1,2,4 and rho_2,5
if Requirement == "M18":
	Req_134, Req_25 = Set_Mandelbaum_Constraints()
elif Requirement == "Z18":
	Req_134 = gaussian_filter( xip_data, Smooth_Scale ) / 10.
	Req_25  = gaussian_filter( xip_data, Smooth_Scale ) / 10.
else:
	print "Requirement must be set to M18 or Z18"
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
					delta_xip = T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alpha*(rho[lfv,1,:]+rho[lfv,4,:])
					ax.errorbar(theta, delta_xip, yerr=rho_err[lfv,1,:], color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])
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
			delta_xip = T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alpha*(rho[lfv,1,:]+rho[lfv,4,:])
			err_delta_xip = ( T_ratio**2*(rho_err[lfv,0,:]**2+rho_err[lfv,2,:]**2+rho_err[lfv,3,:]**2) - T_ratio*alpha*(rho_err[lfv,1,:]**2+rho_err[lfv,4,:]**2) )**0.5
			ax.errorbar(theta, abs(delta_xip), yerr=abs(err_delta_xip), color='black', linewidth=3)
			ax.set_xlabel(r'$\theta$ [arcmin]')

	plt.subplots_adjust(hspace=0, wspace=0)
	plt.savefig('LFver%s/rho1/Plot_Overall-rho%s_CovPatches%sx%s_Require%s_%sxip.png'%(LFver[lfv],pm,Res,Res, Requirement, Which_Data))
	#plt.show()
	return
Plot_2Panel(rhop_mean, rhop_err, '+', -1)




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
	delta_xip = T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alpha*(rho[lfv,1,:]+rho[lfv,4,:])
	# append to itself 15 times for the 15 bins (delta_xip same for every bin).
	# NB: THIS NEEDS TO BE CHANGED, TO HAVE ALPHA VALUES AND THUS DELTA_XIP SPECIFIC TO EACH BIN.
	delta_xip_stack = np.tile(delta_xip, 15)

	# Assemble the theory vector - used to guage signif. of measuring genuine signal.
	# And deviations in this signal from those with high/low values of S_8
	def Read_In_Theory_Vector(hi_lo_fid):
		# hi_lo_fid must be one of 'high', 'low' or 'fid'
		indir_theory = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/ForBG/outputs/test_output_S8_%s_test/shear_xi_plus/' %hi_lo_fid
		theta_theory = np.loadtxt('%s/theta.txt' %indir_theory) * (180./np.pi) * 60.	# Convert long theta array in radians to arcmin
		xip_theory_stack = np.zeros( [15,len(theta)] )									# Will store all auto & cross xi_p for the 5 tomo bins
		idx = 0
		for i in range(1,6):
			for j in range(1,6):
				if i >= j:		# Only read in bins 1-1, 2-1, 2-2, 3-1, 3-2,...
					tmp_xip_theory = np.loadtxt('%s/bin_%s_%s.txt' %(indir_theory,i,j))
					xip_theory_stack[idx,:] = np.interp( theta, theta_theory, tmp_xip_theory )		# sample at the theta values used for PSF modelling.
					idx+=1
		xip_theory_stack = xip_theory_stack.flatten()
		return xip_theory_stack
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
		chi2_sys[i] = np.dot( np.transpose(noise+delta_xip_stack), 
								np.dot( np.linalg.inv(Cov_Mat_uc_Survey_All), noise+delta_xip_stack ))
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
print " It took %.0f seconds." %(t2-t1)



