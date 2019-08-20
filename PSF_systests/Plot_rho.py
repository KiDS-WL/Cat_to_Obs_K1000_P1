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
data_dir = '/disk2/ps1/bengib/KiDS1000_NullTests/Codes_4_My_Eyes/Catalogues/'
data_ZBlabel = '0.1-0.3'		# '0.1-0.3' 'None'
# North
theta_data, xip_data_N_tmp, npairs_data_N = np.loadtxt('%s/KN.BlindA.xi_pm.ZBcut%s.dat'%(data_dir,data_ZBlabel), usecols=(0,1,3), unpack=True)
xip_data_N = np.interp(theta, theta_data, xip_data_N_tmp)
# South
theta_data, xip_data_S_tmp, npairs_data_S = np.loadtxt('%s/KS.BlindA.xi_pm.ZBcut%s.dat'%(data_dir,data_ZBlabel), usecols=(0,1,3), unpack=True)
xip_data_S = np.interp(theta, theta_data, xip_data_S_tmp)
xip_data = (xip_data_N*npairs_data_N + xip_data_S*npairs_data_S) / (npairs_data_N + npairs_data_S)
if data_ZBlabel == '0.1-0.3':
	# Lowest Zb bin has a wee bit of noise in it. 
	# Smoothing it a wee bit for plotting purposes.
	sigma = 0.8
	xip_data = gaussian_filter( xip_data, sigma )


def Set_Scales(ax):
	ax.set_xscale('log')
	ax.set_xlim([0.5,300.])
	#ax.set_yscale('log')
	return

# LIMITS = [ [-13.9,2], [0.1,40], [-0.04,0.25], [-0.29,0.1], [-25,-0.0], [-15, 0.9] ]
LIMITS = [ [-13.9,40], [-13.9,40], [-0.29,0.25], [-0.29,0.25], [-25, 1.9], [-25, 1.9] ]

def Plot_5(rho, rho_err, pm):

	Deltaxip = True
	if Deltaxip:
		scrollthrough=6
	else:
		scrollthrough=5

	# Define quantities needed for delta_xip
	alpha=0.03 		# worse case scenario PSF leakage is 0.03;
			   		# but you still need to calculate this term from data
	T_ratio = 1. 	# Will ultimately get this column from Lance's updated Lensfit cats [I think?]

	fig = plt.figure(figsize = (24,10))
	gs1 = gridspec.GridSpec(3, 2)
	colors = [ 'magenta', 'darkblue', 'dimgrey', 'orange', 'lawngreen', 'cyan' ] 	
	for i in range(scrollthrough):
			ax = plt.subplot(gs1[i], adjustable='box')
			Set_Scales(ax)
			ax.set_ylim( LIMITS[i] )
		
			ax.fill_between(theta_data[:], y1=-0.1*xip_data*1e7, y2=0.1*xip_data*1e7, facecolor='yellow') 
			if i ==5:
				ax.set_ylabel(r'$\delta\xi_%s(\theta) \times 10^{-7}$'%pm)
				ax.set_xlabel(r'$\theta$ [arcmin]')
				for lfv in range(len(LFver)):
					delta_xip = T_ratio**2*(rho[lfv,0,:]+rho[lfv,2,:]+rho[lfv,3,:]) - T_ratio*alpha*(rho[lfv,1,:]+rho[lfv,4,:])
					ax.errorbar(theta, delta_xip*1e7, yerr=rho_err[lfv,1,:]*1e7, color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])
			else:
				ax.set_ylabel(r'$\rho_{%s}(\theta) \times 10^{-7}$'%(i+1))
				if i!=3 and i!=4:
					ax.set_xticks([])
				else:
					ax.set_xlabel(r'$\theta$ [arcmin]')
				for lfv in range(len(LFver)):
					ax.errorbar(theta, rho[lfv,i,:]*1e7, yerr=rho_err[lfv,i,:]*1e7, color=colors[lfv], linewidth=3, label=r'%s'%Plot_Labels[lfv])


	ax.legend(loc='lower right', frameon=False)
	#fig.suptitle('Lensfit v'+LFver[0])
	if Deltaxip == False:
		plt.subplot(gs1[-1], adjustable='box').set_visible(False)
	plt.subplots_adjust(hspace=0)
	plt.savefig('LFver%s/rho1/Plot_rho%s_CovPatches%sx%s.png'%(LFver[0],pm,Res,Res))
	#plt.show()
	return

#Plot_5(rhop_mean, rhop_err, '+')
#Plot_5(rhom_mean, rhom_err, '-')


def Plot_4Paper(rho, rho_err, pm, lfv):
	# Define quantities needed for delta_xip
	alpha=0.03 		# worse case scenario PSF leakage is 0.03;
			   		# but you still need to calculate this term from data
	T_ratio = 1. 	# Will ultimately get this column from Lance's updated Lensfit cats.

	fig = plt.figure(figsize = (10,10))
	gs1 = gridspec.GridSpec(2, 1)
	colors = [ 'magenta', 'darkblue', 'dimgrey', 'orange', 'lawngreen', 'cyan' ] 	
	for i in range(2):
		ax = plt.subplot(gs1[i], adjustable='box')
		Set_Scales(ax)
		ax.set_yscale('log')
		ax.fill_between(theta_data[:], y1=-0.1*xip_data, y2=0.1*xip_data, facecolor='yellow') 
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
	plt.savefig('LFver%s/rho1/Plot_Overall-rho%s_CovPatches%sx%s.png'%(LFver[lfv],pm,Res,Res))
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
	alpha=0.03 	# worse case scenario PSF leakage is 0.03;
			   	# but you still need to calculate this term from data
	T_ratio = 1. # Will ultimately get this column from Lance's updated Lensfit cats.
	
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




