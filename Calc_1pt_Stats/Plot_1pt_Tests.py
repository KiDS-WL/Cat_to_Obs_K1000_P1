# 18/10/2018, B. M. Giblin, PhD student Edinburgh
# Just plot some stuff with the K1000 preliminary catalogues
# CH: 10th Sept - update to metacal 'mc' and autocal 'ac' cats

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from scipy.stats.mstats import mquantiles
from scipy.stats import binned_statistic_2d
import sys
import os 
import time

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}

plt.rc('font', **font)

# Read in user input to set the patch, blind, zmin,zmax, nbootstrap
if len(sys.argv) <5: 
    print "Usage: %s Patch Blind ZBmin ZBmax nbootstrap" % sys.argv[0] 
    sys.exit(1)
else:
	NorS=sys.argv[1]  # N for North, S for South, or 'All' for both.
	Blind=sys.argv[2] # Just put 'A' 	
	ZBlo=np.float(sys.argv[3])  # If numbers are inputted, this Z_B cut is applied to data
	ZBhi=np.float(sys.argv[4])  # For no ZBcut, put 'None' for both of these
	nboot=np.int(sys.argv[5]) # this is the number of bootstrap realisations to run


Read_Cat_Or_Pickle = "Cat" 
#Read_Cat_Or_Pickle = "Pickle" 	# Set to "Cat" to read from fitsfile catalogue (~680s for N&S)
								# Or "Pickle" to read smaller pickled version (~20s for N&S)

DIRECT = os.getcwd() 

def Pickle_Data(ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR):
	print "Pickling RA,Dec,e1,e2,w,ZB...."
	outputrd='%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB_all' %(DIRECT,NorS,Blind)
	np.save(outputrd, np.column_stack((ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB)) )

	print "Pickling Xpos,Ypos,MAG,PSFe1,PSFe2,TPSF,SNR,ZB...."
	outputXY='%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB_all' %(DIRECT,NorS,Blind)
	np.save(outputXY, np.column_stack((Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR,ZB)) )

	return


def Read_Basic_Data(Read_NorS):
	print "Reading in RA,Dec,e1,e2,w,ZB...."

	# This is currently v messy - sorry Ben - we will clean this up when we have a stable catalogue and format
	# KV450 cats
	f = fits.open('/disk10/cech/KIDSCOLLAB/KIDSCOLLAB_V0.5.9/KV450_CATALOGUES_PATCH_V0.5.9/KV450_%s_reweight_3x4x4_v2_good.cat'%(Read_NorS))
	# fits extension for KV450 G9, G12, G15 etc iext=1, for GAll
	iext=1

	# K1000 mc cats
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/K1000_%s_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321.cat'%(Read_NorS))
	#fits extension for K1000
	#iext=2

	# K1000-mc preeblind cats
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/recal_split/K1000_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_%s.cat'%(Read_NorS))
	#fits extension for K1000
	#iext=2

	# K1000- Phase 0 cats
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/K1000_%s_9band_mask_BLINDED_v3.cat'%(Read_NorS))
	#fits extension for K1000
	#iext=2

	# K1000- G9 Frankenstein cats
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/%s_9band_mask_LF_svn_309b_FRANKENSTEIN.cat'%(Read_NorS))
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/%s_9band_mask_LF_svn_309b.cat'%(Read_NorS))
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/%s_9band_mask_LF_glab_321.cat'%(Read_NorS))
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/%s_9band_mask_LF_glab_321_FRANKENSTEIN.cat'%(Read_NorS))
	#f = fits.open('/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/%s_9band_mask_LF_glab_322.cat'%(Read_NorS))
	#fits extension for K1000
	#iext=2

	ra = f[iext].data['ALPHA_J2000']
	dec = f[iext].data['DELTA_J2000']
	ZB = f[iext].data['Z_B']

	
	#KV450 columns - paste autocal into metacal to use K1000 script
	mce1 = f[iext].data['bias_corrected_e1']
	mce2 = f[iext].data['bias_corrected_e2']
	#mcw = f[iext].data['recal_weight']
	mcw = f[iext].data['weight']
	ace1 = f[iext].data['bias_corrected_e1'] - f[iext].data['e1_correction']
	ace2 = f[iext].data['bias_corrected_e2'] - f[iext].data['e2_correction']
	#acw = f[iext].data['recal_weight']
	acw = f[iext].data['weight']
	ngals = len(acw)
	acm1 = np.ones(ngals)
	acm2 = np.ones(ngals)
	mcm1 = np.ones(ngals)
	mcm2 = np.ones(ngals)
	Xpos = f[iext].data['Xpos_THELI']
	Ypos = f[iext].data['Ypos_THELI']
    
	'''
	#K1000 columns Phase 0
	#metacal columns
	mce1 = f[iext].data['e1_%s'%Blind]
	mce2 = f[iext].data['e2_%s'%Blind]
	mcw = f[iext].data['weight_%s'%Blind]
	#autocal columns
	ace1 = f[iext].data['e1_%s'%Blind]
	ace2 = f[iext].data['e2_%s'%Blind]
	acw = f[iext].data['weight_%s'%Blind]
	ngals = len(acw)
	acm1 = np.ones(ngals)
	acm2 = np.ones(ngals)
	mcm1 = np.ones(ngals)
	mcm2 = np.ones(ngals)
	'''
	'''
	#K1000 columns
	#metacal columns
	mce1 = f[iext].data['raw_e1_%s'%Blind]
	mce2 = f[iext].data['raw_e2_%s'%Blind]
	mcw = f[iext].data['metacal_weight_%s'%Blind]
	mcm1 = 1+f[iext].data['metacal_m1']
	mcm2 = 1+f[iext].data['metacal_m2']
	#autocal columns
	ace1 = f[iext].data['autocal_e1_%s'%Blind]
	ace2 = f[iext].data['autocal_e2_%s'%Blind]
	acw = f[iext].data['recal_weight_%s'%Blind]
	ngals = len(acw)
	acm1 = np.ones(ngals)
	acm2 = np.ones(ngals)
	'''
	'''
	#K1000 columns - preblind, pre mc weight calc
	#metacal columns
	mce1 = f[iext].data['raw_e1']
	mce2 = f[iext].data['raw_e2']
	mcm1 = 1+f[iext].data['metacal_m1']
	mcm2 = 1+f[iext].data['metacal_m2']
	mcw=(mcm1+mcm2)*0.5
	#autocal columns
	ace1 = f[iext].data['autocal_e1']
	ace2 = f[iext].data['autocal_e2']
	#acw = f[iext].data['recal_weight']
	acw = f[iext].data['weight']
	ngals = len(acw)
	acm1 = np.ones(ngals)
	acm2 = np.ones(ngals)
	acw=np.ones(ngals)
	'''

	print "Reading in Xpos,Ypos,MAG,PSFe's...."
	#Xpos = f[iext].data['Xpos']
	#Ypos = f[iext].data['Ypos']
	MAG = f[iext].data['MAG_AUTO']
	PSFe1 = f[iext].data['PSF_e1']
	PSFe2 = f[iext].data['PSF_e2']
	Ixx = f[iext].data['PSF_Q11']
	Iyy = f[iext].data['PSF_Q22']
	TPSF = f[iext].data['PSF_Q11']*f[iext].data['PSF_Q22'] - f[iext].data['PSF_Q12']*f[iext].data['PSF_Q12']
	SNR = f[iext].data['pixel_SNratio']
	f.close()
	return ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR


def Stack_NandS():
	print "Stacking North and South..."
	t1 = time.time()
	ra = np.append(ra_N, ra_S)
	dec = np.append(dec_N, dec_S)
	mce1 = np.append(mce1_N, mce1_S)
	mce2 = np.append(mce2_N, mce2_S)
	mcm1 = np.append(mcm1_N, mcm1_S)
	mcm2 = np.append(mcm2_N, mcm2_S)
	mcw = np.append(mcw_N, mcw_S)
	ace1 = np.append(ace1_N, ace1_S)
	ace2 = np.append(ace2_N, ace2_S)
	acw = np.append(acw_N, acw_S)
	ngals=len(acw)
	acm1=np.ones(ngals)
	acm2=np.ones(ngals)
	ZB = np.append(ZB_N, ZB_S)

	Xpos = np.append(Xpos_N, Xpos_S)
	Ypos = np.append(Ypos_N, Ypos_S)
	MAG = np.append(MAG_N, MAG_S)
	PSFe1 = np.append(PSFe1_N, PSFe1_S)
	PSFe2 = np.append(PSFe2_N, PSFe2_S)
	TPSF = np.append(TPSF_N, TPSF_S)
	SNR = np.append(SNR_N, SNR_S)
	t2 = time.time()
	return ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR


if Read_Cat_Or_Pickle == "Cat":


	t1 = time.time()
	if NorS=="All":
		ra_N,dec_N,mce1_N,mce2_N,mcw_N,mcm1_N,mcm2_N,ace1_N,ace2_N,acw_N,ZB_N, Xpos_N,Ypos_N,MAG_N, PSFe1_N,PSFe2_N,TPSF_N,SNR_N = Read_Basic_Data("N")
		ra_S,dec_S,mce1_S,mce2_S,mcw_S,mcm1_S,mcm2_S,ace1_S,ace2_S,acw_S,ZB_S, Xpos_S,Ypos_S,MAG_S, PSFe1_S,PSFe2_S,TPSF_S,SNR_S = Read_Basic_Data("S")
		ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Stack_NandS()
	else:
		ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Read_Basic_Data(NorS)
	t2 = time.time()
	print "It took %.0f seconds to read in from the catalogue" %(t2-t1)	
	Pickle_Data(ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR)
	
	# set the autocal m-calibration correction to 1 - not included in the pickle files for speed
	ngals=len(ra)				
	acm1=np.ones(ngals)
	acm2=np.ones(ngals)	



elif Read_Cat_Or_Pickle == "Pickle":


	t1 = time.time()
	if NorS=="All":

		# Check if pickled 'All' catalogue exists... If so, proceed...
		if os.path.isfile('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB_all.npy'%(DIRECT,NorS,Blind)):
			print "Reading in RA,Dec,e1,e2,w,ZB...."
			ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB = np.load('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB_all.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6,7,8,9,10],:]
			print "Reading in Xpos,Ypos,MAG,PSFe's...."
			Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = np.load('%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB_all.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6],:]

		# ...if pickled 'All; catalogue does not exist, make it, to save time stacking North and South next time...
		else:
			print "Reading in RA,Dec,e1,e2,w,ZB...."
			ra_N,dec_N,mce1_N,mce2_N,mcw_N,mcm1_N,mcm2_N,ace1_N,ace2_N,acw_N,ZB_N = np.load('%s/Catalogues/KN.Blind%s.ra_dec_e1_e2_w_ZB_all.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6,7,8,9,10],:]
			ra_S,dec_S,mce1_S,mce2_S,mcw_S,mcm1_S,mcm2_S,ace1_S,ace2_S,acw_S,ZB_S = np.load('%s/Catalogues/KS.Blind%s.ra_dec_e1_e2_w_ZB_all.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6,7,8,9,10],:]
			print "Reading in Xpos,Ypos,MAG,PSFe's...."
			Xpos_N,Ypos_N,MAG_N,PSFe1_N,PSFe2_N,TPSF_N,SNR_N = np.load('%s/Catalogues/KN.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB_all.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6],:]
			Xpos_S,Ypos_S,MAG_S,PSFe1_S,PSFe2_S,TPSF_S,SNR_S = np.load('%s/Catalogues/KS.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB_all.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6],:]
			ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Stack_NandS()
			Pickle_Data(ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR)

	else:
		print "Reading in RA,Dec,e1,e2,w,ZB...."
		ra,dec,mce1,mce2,mcw,mcm1,mcm2,ace1,ace2,acw,ZB = np.load('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB_all.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6,7,8,9,10],:]
		print "Reading in Xpos,Ypos,MAG,PSFe's...."
		Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = np.load('%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB_all.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6],:]

		

	t2 = time.time()
	print "It took %.0f seconds to read in from the pickled files" %(t2-t1)
	
	# set the autocal m-calibration correction to 1 - not included in the pickle files for speed
	ngals=len(ra)
	acm1=np.ones(ngals)
	acm2=np.ones(ngals)	

	print np.min(TPSF), np.max(TPSF)

if str(ZBlo) == "None":
    ZBlabel = 'ZBcutNone'
elif float(ZBlo)<float(ZBhi):
	print "Making the ZBcut in the range %s to %s" %(ZBlo, ZBhi)
	ZBlabel = 'ZBcut%s-%s' %(ZBlo,ZBhi)
	#idx = np.where( np.logical_and(ZB>float(ZBlo), ZB<float(ZBhi)) )[0]
	idx=( (ZB>ZBlo) & (ZB<ZBhi) & (TPSF>4))
	ra=ra[idx]
	dec=dec[idx]
	#metacal columns
	mce1=mce1[idx]
	mce2=mce2[idx]
	mcw=mcw[idx]
	mcm1=mcm1[idx]
	mcm2=mcm2[idx]
	#autocal columns
	ace1=ace1[idx]
	ace2=ace2[idx]
	acw=acw[idx]
	acm1=acm1[idx]
	acm2=acm2[idx]
	#Other
	ZB=ZB[idx]
	Xpos=Xpos[idx]
	Ypos=Ypos[idx]
	MAG=MAG[idx]
	PSFe1=PSFe1[idx]
	PSFe2=PSFe2[idx]
	TPSF=TPSF[idx]
	SNR=SNR[idx]
else:
	print "ZBlo %s is not lower than ZBhi %s. Amend this or set them both to None."
	sys.exit()
# Save the weighted additive bias - and the weighted m-correction to apply to the bias (in the metacal case)
np.savetxt('%s/Catalogues/K%s.Blind%s.ccorr_e1_e2.%s.dat'%(DIRECT,NorS,Blind,ZBlabel), 
           np.c_[ np.average(mce1,weights=mcw),np.average(mce2,weights=mcw),
		          np.average(mcm1,weights=mcw),np.average(mcm2,weights=mcw),
		          np.average(ace1,weights=acw),np.average(ace2,weights=acw) ])


def Bootstrap_Error(nboot, samples, weights,mcorr):
	N = len(samples)
	bt_samples = np.zeros(nboot)		 		# Will store mean of nboot resamples
	for i in range(nboot):
		idx = np.random.randint(0,N,N)			# Picks N random indicies with replacement
		bt_samples[i] = np.sum( weights[idx]*samples[idx] ) / np.sum( weights[idx]*mcorr[idx] )
	return np.std( bt_samples )


def MeanQ_VS_XY(Q, w, m, X,Y,num_XY_bins):
	# we want the weighted mean of Q and also calibrated.  Calculate the sum 2D binned value of Q*w and m*w, and then divide

	sumQw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q*w, statistic='sum', bins=num_XY_bins)#, range=[[0.0,21000.],[0.0,21000.0]])
	sum_mw_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, m*w, statistic='sum', bins=num_XY_bins)#, range=[[0.0,21000.],[0.0,21000.0]])
	
	AvQ_grid=sumQw_grid/sum_mw_grid

	return AvQ_grid,yedges,xedges

def Plot_XY_Grids(Q, yedges,xedges, lowlim, upplim, label, savename):
	cmap = plt.cm.rainbow
	cmap.set_bad('w', 1.)
	plt.figure(figsize=(8,7))
	plt.imshow(Q, cmap=cmap, interpolation='none', vmin=lowlim, vmax=upplim, origin='lower', 
					extent=[np.min(yedges), np.max(yedges),
                            np.min(xedges), np.max(xedges)],aspect='equal')  
	plt.ylabel(r'$Y_{\rm{pos}}$')
	plt.xlabel(r'$X_{\rm{pos}}$')
	plt.colorbar(label=label)
	plt.tight_layout()
	plt.savefig(savename)
	plt.show()
	return

def Plot_BinQx_VS_BinQy(Qx, Qy, weights, mcorr, num_bins, labels, xlabel, ylabel, title,xlimits, ylimits, savename, Bootstrap):

	colors = ['magenta', 'dimgrey', 'darkblue', 'orange', 'lawngreen', 'red']
	f, ((ax1)) = plt.subplots(1, 1, figsize=(8,7))

	if len(Qx.shape) == 1:
		Qx = np.reshape(Qx, (1,len(Qx)))
		Qy = np.reshape(Qy, (1,len(Qy)))

	# Average the corresponding Y's to the X's
	binx_centres = np.zeros([Qy.shape[0], num_bins])
	biny_centres = np.zeros([Qy.shape[0], num_bins])
	biny_centres_err = np.zeros([Qy.shape[0], num_bins])
	probability = np.linspace(0,1,num_bins+1)	

	for i in range(Qy.shape[0]):			
		binx_edges = mquantiles(Qx[i,:], prob=probability)				# Edges of the num_bins bins, each containing same-% of data points.
		
		for j in range(num_bins):
			idx = np.where(np.logical_and(Qx[i,:]>=binx_edges[j], Qx[i,:]<binx_edges[j+1]))[0]
			# print len(idx) / float(len(Qx[i,:]))
			binx_centres[i,j] = np.sum( weights[idx]*Qx[i,idx] ) / np.sum( weights[idx] ) 
			biny_centres[i,j] = np.sum( weights[idx]*Qy[i,idx] ) / np.sum( weights[idx] *mcorr[idx] ) 
			#The number of bootstrap resamples to take in estimating the errors is set on the command line.
			if Bootstrap:
				biny_centres_err[i,j] = Bootstrap_Error(nboot, Qy[i,idx], weights[idx],mcorr[idx])
   			else:
				biny_centres_err[i,j] = 0.
 
		plt.errorbar(binx_centres[i,:], biny_centres[i,:], yerr=biny_centres_err[i,:], fmt='o', color=colors[i],label=labels[i])

		#you might want to add in an offset e.g
		#offset = [-0.01,0.01, -0.03,0.03] 
		# +offset[i]*np.mean(binx_centres)

	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	if ylimits != None:
		plt.ylim(ylimits)
	if xlimits != None:
		plt.plot( xlimits, [0.,0.], 'k--' )
		plt.xlim(xlimits)
	plt.legend(loc='best', frameon=False)
	plt.savefig(savename)
	plt.show()
	return np.dstack(( binx_centres, biny_centres, biny_centres_err ))


# ------------------------------------------------------------- 1-POINT STATISTICS PLOTS ------------------------------------------------------------- #
# NOTE ON RETURNED ARRAY: bin_x_y_yerr[:,:,i] --> i=0 for x, i=1 for y, i=2 for yerr.

#print "Exiting as the next 1pt plots are expensive (~100s) to produce due to bootstrapping errors."
#sys.exit()

Bootstrap = False	# Bootstrap errors? Takes ~nboot times longer! (nboot defined on command line)

for mc_or_ac in ('autocal','metacal'):
	if mc_or_ac == 'metacal':
		e1=mce1
		e2=mce2
		w=mcw
		m=(mcm1+mcm2)*0.5
		m1=mcm1
		m2=mcm2
	else:
		e1=ace1
		e2=ace2
		w=acw
		ngals=len(e1)
		m=np.ones(ngals)
		m1=np.ones(ngals)
		m2=np.ones(ngals)

	'''
	#### e vs X-Y
	# Weighted mean-e in bins VS Xpos,Ypos
	num_XY_bins=20
	e1_grid,yedges,xedges = MeanQ_VS_XY(e1, w, m1, Xpos,Ypos, num_XY_bins)
	e2_grid,yedges,xedges = MeanQ_VS_XY(e2, w, m2, Xpos,Ypos, num_XY_bins)
	#scale the result
	e1_grid=e1_grid*1e3
	e2_grid=e2_grid*1e3
	#auto-set the cmap scale to span 80% of the min/max range
	#v1min = np.min(e1_grid) + 0.1*(np.max(e1_grid) - np.min(e1_grid))
	#v1max = np.max(e1_grid) - 0.1*(np.max(e1_grid) - np.min(e1_grid))
	#v2min = np.min(e2_grid) + 0.1*(np.max(e2_grid) - np.min(e2_grid))
	#v2max = np.max(e2_grid) - 0.1*(np.max(e2_grid) - np.min(e2_grid))
	#or maybe you want a fixed range to compare easily
	v1min=-8.5
	v1max=9.4
	v2min=v1min
	v2max=v1max

	Plot_XY_Grids(e1_grid, yedges,xedges,v1min, v1max, r'$\langle e_1 \rangle\times 10^{-3}$', '%s/GeneralPlots/K%s.%s.Blind%s.Meane1_VS_XY.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel) )
	Plot_XY_Grids(e2_grid, yedges,xedges,v2min, v2max, r'$\langle e_2 \rangle\times 10^{-3}$', '%s/GeneralPlots/K%s.%s.Blind%s.Meane2_VS_XY.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel) )
	'''

	# MEAN ELLIPTICITY VS X/Ypos
	#bin_Xpos_e_eerr  = Plot_BinQx_VS_BinQy(np.vstack((Xpos,Xpos)), np.vstack((e1, e2)), w, 10, 
	#[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$X_{\rm{pos}}$', r'$\langle e \rangle$', 
	#None, '%s/GeneralPlots/K%s.Blind%s.Xpos_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )

	#bin_Ypos_e_eerr = Plot_BinQx_VS_BinQy(np.vstack((Ypos,Ypos)), np.vstack((e1, e2)), w, 10, 
	#[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$Y_{\rm{pos}}$', r'$\langle e \rangle$', 
	#None, '%s/GeneralPlots/K%s.Blind%s.Ypos_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )

	t1 = time.time()
	# MEAN ELLIPTICITY VS PSF ELLIPTICITY
	bin_PSFe_e_eerr = Plot_BinQx_VS_BinQy(np.vstack((PSFe1,PSFe2)), 1000*np.vstack((e1, e2)), w, m, 15, 
				[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$e_{\rm{PSF}}$', r'$\langle e \rangle \times 10^{-3}$', mc_or_ac,
				[-0.04,0.04], [-7.5,7.5], 
				'%s/GeneralPlots/K%s.%s.Blind%s.ePSF_VS_e.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel), Bootstrap )
	t2 = time.time()
	print "It took %.0f s to bin and plot e VS PSFe" %(t2-t1)
    
	# MEAN ELLIPTICITY VS ZB
	#bin_ZB_e_err = Plot_BinQx_VS_BinQy(np.vstack((ZB,ZB)), 1000*np.vstack((e1, e2)), w, m, 12, 
	#				[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$z_{\rm B}$', r'$\langle e \rangle \times 10^{-3}$', mc_or_ac,
	#			[0.1,2.0], [-2.5,3.5], 
	#			'%s/GeneralPlots/K%s.%s.Blind%s.ZB_VS_e.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel), Bootstrap  )
	#print "Finished binning and plotting e VS ZB" 

	# MEAN ELLIPTICITY VS MAGNITUDE
	#bin_MAG_e_err = Plot_BinQx_VS_BinQy(np.vstack((MAG,MAG)), 1000*np.vstack((e1, e2)), w, m,12, 
	#				[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'Magnitude', r'$\langle e \rangle \times 10^{-3}$',mc_or_ac,
	#				[20,25.5], [-5.0,5.0], 
	#			'%s/GeneralPlots/K%s.%s.Blind%s.MAG_VS_e.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel), Bootstrap  )
	#print "Finished binning and plotting e VS MAG" 


	# MEAN ELLIPTICITY VS PSF SIZE
	#bin_TPSF_e_err = Plot_BinQx_VS_BinQy(np.vstack((TPSF,TPSF)), 1000*np.vstack((e1, e2)), w, m,12, 
	#				[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$T_{\rm{PSF}}$', r'$\langle e \rangle \times 10^{-3}$',mc_or_ac,	
	#				[2.5,5.5], [-2.5,2.5], 
	#				'%s/GeneralPlots/K%s.%s.Blind%s.TPSF_VS_e.%s.png' %(DIRECT,NorS,mc_or_ac,Blind,ZBlabel), Bootstrap  )
	#print "Finished binning and plotting e VS TPSF" 





'''
# May want to keep this general function for making a histogram of the data.

def Histogram(Q, num_bins, labels, xlabel, title, savename):

    colors = ['magenta', 'dimgrey', 'darkblue', 'orange', 'lawngreen', 'red']
    f, ((ax1)) = plt.subplots(1, 1, figsize=(9,9))
    if len(Q.shape) == 1:
        Q = np.reshape(Q, (1,len(Q)))
        
    for i in range(Q.shape[0]):
        histo, tmp_bins = np.histogram(Q[i,:], num_bins)
        # get bin centres
        bins = tmp_bins[:-1] + (tmp_bins[1] - tmp_bins[0])/2.
        ax1.bar(bins, histo, width=(bins[1]-bins[0]), color=colors[i], edgecolor=colors[i], alpha=0.5, label=labels[i])
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(r'PDF')
    ax1.legend(loc='best', frameon=False)
    ax1.set_title(title)
    plt.savefig(savename)
    plt.show()
    return

#Histogram(np.vstack((e1,e2)), 100, [r'$e_1$',r'$e_2$'], r'$e$', r'KiDS-%s'%NorS, '%s/GeneralPlots/K%s.Blind%s.eHist.%s.png' %(DIRECT,NorS,Blind,ZBlabel) )
'''







