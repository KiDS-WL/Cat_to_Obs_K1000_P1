import numpy as np
import pylab as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter
import sys
import os 
import glob
rc('text',usetex=True)
rc('font',size=24)
rc('legend',**{'fontsize':24})
rc('font',**{'family':'serif','serif':['Computer Modern']})

# Set the plotted requirement on the rho stats:
Requirement = "M18"			# "M18" for Mandelbaum+18: [p_1,2,4<xi_+/SNR, p_2,5<xi_+/(SNR*alpha)
							# "Z18" for Zuntz+18: xi_+ / 10 

Which_Data = "Data"		# "Data" for K1000 or "Theory" for HaloFit
							# goes into calculation of the requirement bands.

# Define quantities needed for delta_xip
alpha=0.03 		# worse case scenario PSF leakage is 0.03;
			   	# but you still need to calculate this term from data
T_ratio = 1. 	# Will ultimately get this column from Lance's updated Lensfit cats [I think?]

# If using the lowest ZB bin, you might want to smooth the data vector 
# so the Requirement band is less jaggardy. If so, use this smoothign scale:
Smooth_Scale = 0.8


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
rhom_mean = np.zeros_like(rhop_mean)
# Get error (on mean) bars
rhop_err = np.zeros_like(rhop_mean)
rhom_err = np.zeros_like(rhom_mean)

for lfv in range(len(LFver)):

	rhop_split = np.zeros([ numN[lfv]+numS[lfv], 5, ThBins ]) 
	rhom_split = np.zeros_like(rhop_split)
	#meanfactor = 2. # Divide by this, unless no N or S Field, in which case set to 1.
	for i in range(5):
		try:
			theta, rhopN, rhomN = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_N.dat'%(LFver[lfv],i+1,i+1), usecols=(0,1,2), unpack=True)
			# If the above exists, try to read in the weight (only saved this for LFver321)
			try:
				weightN = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_N.dat'%(LFver[lfv],i+1,i+1), usecols=(3,), unpack=True)
			except IndexError:
				weightN = 1.
		except IOError:
			weightN = 1.
			rhopN = 0.
			rhomN = 0.
		try:
			theta, rhopS, rhomS = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_S.dat'%(LFver[lfv],i+1,i+1), usecols=(0,1,2), unpack=True)
			# If the above exists, try to read in the weight (only saved this for LFver321)
			try:
				weightS = np.loadtxt('LFver%s/rho%s/rho%s_KiDS_S.dat'%(LFver[lfv],i+1,i+1), usecols=(3,), unpack=True)
			except IndexError:
				weightS = 1.
		except IOError:
			weightS = 1.
			rhopS = 0.
			rhomS = 0.
	
		# Weighted average of rho-North & rho-South
		rhop_mean[lfv,i,:] = (weightN*rhopN + weightS*rhopS) / (weightN+weightS)
		rhom_mean[lfv,i,:] = (weightN*rhomN + weightS*rhomS) / (weightN+weightS)
		for j in range(numN[lfv]):
			rhop_split[j,i,:], rhom_split[j,i,:] = np.loadtxt(NFiles[lfv][j], usecols=(1,2), unpack=True)
		for j in range(numS[lfv]):
			rhop_split[numN[lfv]+j,i,:], rhom_split[numN[lfv]+j,i,:] = np.loadtxt(SFiles[lfv][j], usecols=(1,2), unpack=True)

		rhop_err[lfv,i,:] = np.sqrt( np.diag( np.cov(rhop_split[:,i,:], rowvar = False) ) / (numN[lfv]+numS[lfv]) ) 
		rhom_err[lfv,i,:] = np.sqrt( np.diag( np.cov(rhom_split[:,i,:], rowvar = False) ) / (numN[lfv]+numS[lfv]) ) 

# To get a rough idea of size of rho stats, read in the xi+- of some data to overplot 
data_dir = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/xi_pm/'
data_ZBlabel = '0.1-0.3'		# '0.1-0.3' 'None'
theta_data, xip_data = np.loadtxt('%s/KAll.BlindA.xi_pm.ZBcut%s.dat' %(data_dir, data_ZBlabel), usecols=(0,1), unpack=True)
 
theta_theory, xip_theory_tmp = np.loadtxt('%s/xi_p_smith03revised_zKV450_ZBcut%s' %(data_dir, data_ZBlabel), usecols=(0,1), unpack=True)
xip_theory = np.interp( theta_data, theta_theory, xip_theory_tmp )


def Set_Mandelbaum_Constraints():
	NLOS_Cov = 1250
	cosmol_Cov = 'fid'
	SurveySize = 1000.

	# NOTE BIG DIFFERENCES IN SNR WITH REDSHIFT CUT
	# Linc Cov: 21 (0.1-0.3), 76 (0.7-0.9)
	Cov_inDIR = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/Lincs_CovMat/'
	Cov_Mat_uc_Survey = np.loadtxt('%s/Raw_Cov_Mat_Values.dat' %Cov_inDIR)[0:9, 0:9]			# [0:9, 0:9] This extracts xi+ Cov in lowest bin

	scale = np.where(theta_data < 72)[0]
	SNR = np.dot( np.transpose(xip_data[scale]), 
			np.dot( np.linalg.inv(Cov_Mat_uc_Survey[scale[0]:scale[-1]+1,scale[0]:scale[-1]+1]), xip_data[scale] ))		

	# ! EITHER USE A THEORY VECTOR FOR REQUIREMENTS, OR USE SMOOTHED DATA VECTOR !
	if Which_Data == 'Data':
		# Lowest Z_B bin has a wee bit of noise in it. 
		# Smoothing it a wee bit for plotting purposes.
		Sm_xip = gaussian_filter( xip_data, Smooth_Scale )
	elif Which_Data == 'Theory':
		Sm_xip = np.copy( xip_theory )

	# IF YOU USE 0.1-0.3 DATA/COV THE REQ IS LARGER ON SMALL SCALES AND SMALLER ON LARGE SCALES RELATIVE TO IF WE USED (0.7-0.9)
	# HERE'S THE RATIO: array([ 3.74749516,  2.30591263,  3.22123142,  2.4665902 ,  0.55896293, -0.41969508,  0.39032223,  0.07134366,  0.79189828])
	# SO IN SHORT. NOT OBVIOUS WHICH ONE ZB TO USE. SO GO WITH EITHER. 
	Requirement_134 = T_ratio**(-2) * Sm_xip / (2*SNR)
	Requirement_25  = T_ratio**(-1) * Sm_xip / (2*SNR*alpha)
	return Requirement_134, Requirement_25, Cov_Mat_uc_Survey


# Set the requirements on rho_1,2,4 and rho_2,5
if Requirement == "M18":
	Req_134, Req_25, Cov_Mat_uc_Survey = Set_Mandelbaum_Constraints()
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


def Plot_5_Symlog(rho, rho_err, pm):

	dummy, dummy, Cov_Mat_uc_Survey = Set_Mandelbaum_Constraints()

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
Plot_5_Symlog(rhop_mean, rhop_err, '+')

def Plot_4Paper(rho, rho_err, pm, lfv):

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
	plt.show()
	return
Plot_4Paper(rhop_mean, rhop_err, '+', -1)





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
	plt.show()

	return
#Plot_xip_Plus_rho(rhop_mean, rhop_err, '+')




