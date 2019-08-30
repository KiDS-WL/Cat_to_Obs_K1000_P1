# 18/10/2018, B. M. Giblin, PhD student Edinburgh
# Just plot some stuff with the K1000 preliminary catalogues

import numpy as np
import pylab as plt
from astropy.io import fits
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.stats.mstats import mquantiles
from scipy.stats import binned_statistic_2d
import sys
import os 
import time
rc('text',usetex=True)
rc('font',size=24)
rc('legend',**{'fontsize':30})
rc('font',**{'family':'serif','serif':['Computer Modern']})

NorS=sys.argv[1]  # N for North, S for South, or 'All' for both.
Blind=sys.argv[2] # Just put 'A' 
ZBlo=sys.argv[3]  # If numbers are inputted, this Z_B cut is applied to data
ZBhi=sys.argv[4]  # For no ZBcut, put None for both of these

Read_Cat_Or_Pickle = "Pickle" 	# Set to "Cat" to read from fitsfile catalogue (~680s for N&S)
								# Or "Pickle" to read smaller pickled version (~20s for N&S)


DIRECT = os.getcwd() 


def Pickle_Data(ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR):
	print "Pickling RA,Dec,e1,e2,w,ZB...."
	outputrd='%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone' %(DIRECT,NorS,Blind)
	np.save(outputrd, np.column_stack((ra,dec,e1,e2,w,ZB)) )

	print "Pickling Xpos,Ypos,MAG,PSFe1,PSFe2,TPSF,SNR,ZB...."
	outputXY='%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB.ZBcutNone' %(DIRECT,NorS,Blind)
	np.save(outputXY, np.column_stack((Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR,ZB)) )
	sys.exit()
	return


def Read_Basic_Data(Read_NorS):
	print "Reading in RA,Dec,e1,e2,w,ZB...."
	f = fits.open('%s/Catalogues/K1000_%s_9band_mask_BLINDED_v3.cat'%(DIRECT,Read_NorS))
	ra = f[2].data['ALPHA_J2000']
	dec = f[2].data['DELTA_J2000']
	e1 = f[2].data['e1_%s'%Blind]
	e2 = f[2].data['e2_%s'%Blind]
	w = f[2].data['weight_%s'%Blind]
	ZB = f[2].data['Z_B']
	
	print "Reading in Xpos,Ypos,MAG,PSFe's...."
	Xpos = f[2].data['Xpos']
	Ypos = f[2].data['Ypos']
	MAG = f[2].data['MAG_AUTO']
	PSFe1 = f[2].data['PSF_e1']
	PSFe2 = f[2].data['PSF_e2']
	Ixx = f[2].data['PSF_Q11']
	Iyy = f[2].data['PSF_Q22']
	TPSF = Ixx + Iyy
	SNR = f[2].data['pixel_SNratio']
	f.close()
	return ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR

def Stack_NandS():
	print "Stacking North and South..."
	t1 = time.time()
	ra = np.append(ra_N, ra_S)
	dec = np.append(dec_N, dec_S)
	e1 = np.append(e1_N, e1_S)
	e2 = np.append(e2_N, e2_S)
	w = np.append(w_N, w_S)
	ZB = np.append(ZB_N, ZB_S)

	Xpos = np.append(Xpos_N, Xpos_S)
	Ypos = np.append(Ypos_N, Ypos_S)
	MAG = np.append(MAG_N, MAG_S)
	PSFe1 = np.append(PSFe1_N, PSFe1_S)
	PSFe2 = np.append(PSFe2_N, PSFe2_S)
	TPSF = np.append(TPSF_N, TPSF_S)
	SNR = np.append(SNR_N, SNR_S)
	t2 = time.time()
	return ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR




if Read_Cat_Or_Pickle == "Cat":


	t1 = time.time()
	if NorS=="N" or NorS=="S":
		ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Read_Basic_Data(NorS)
	elif NorS=="All":
		ra_N,dec_N,e1_N,e2_N,w_N,ZB_N, Xpos_N,Ypos_N,MAG_N, PSFe1_N,PSFe2_N,TPSF_N,SNR_N = Read_Basic_Data("N")
		ra_S,dec_S,e1_S,e2_S,w_S,ZB_S, Xpos_S,Ypos_S,MAG_S, PSFe1_S,PSFe2_S,TPSF_S,SNR_S = Read_Basic_Data("S")
		ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Stack_NandS()
	t2 = time.time()
	print "It took %.0f seconds to read in from the catalogue" %(t2-t1)	
	Pickle_Data(ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR)


elif Read_Cat_Or_Pickle == "Pickle":


	t1 = time.time()
	if NorS=="N" or NorS=="S":
		print "Reading in RA,Dec,e1,e2,w,ZB...."
		ra,dec,e1,e2,w,ZB = np.load('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5],:]
		print "Reading in Xpos,Ypos,MAG,PSFe's...."
		Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = np.load('%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB.ZBcutNone.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6],:]


	elif NorS=="All":

		# Check if pickled 'All' catalogue exists... If so, proceed...
		if os.path.isfile('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy'%(DIRECT,NorS,Blind)):
			print "Reading in RA,Dec,e1,e2,w,ZB...."
			ra,dec,e1,e2,w,ZB = np.load('%s/Catalogues/K%s.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5],:]
			print "Reading in Xpos,Ypos,MAG,PSFe's...."
			Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = np.load('%s/Catalogues/K%s.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB.ZBcutNone.npy' %(DIRECT,NorS,Blind)).transpose()[[0,1,2,3,4,5,6],:]

		# ...if pickled 'All; catalogue does not exist, make it, to save time stacking North and South next time...
		else:
			print "Reading in RA,Dec,e1,e2,w,ZB...."
			ra_N,dec_N,e1_N,e2_N,w_N,ZB_N = np.load('%s/Catalogues/KN.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5],:]
			ra_S,dec_S,e1_S,e2_S,w_S,ZB_S = np.load('%s/Catalogues/KS.Blind%s.ra_dec_e1_e2_w_ZB.ZBcutNone.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5],:]
			print "Reading in Xpos,Ypos,MAG,PSFe's...."
			Xpos_N,Ypos_N,MAG_N,PSFe1_N,PSFe2_N,TPSF_N,SNR_N = np.load('%s/Catalogues/KN.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB.ZBcutNone.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6],:]
			Xpos_S,Ypos_S,MAG_S,PSFe1_S,PSFe2_S,TPSF_S,SNR_S = np.load('%s/Catalogues/KS.Blind%s.Xpos_Ypos_MAG_PSFe1_PSFe2_TPSF_SNR_ZB.ZBcutNone.npy' %(DIRECT,Blind)).transpose()[[0,1,2,3,4,5,6],:]
			ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR = Stack_NandS()
			Pickle_Data(ra,dec,e1,e2,w,ZB, Xpos,Ypos,MAG, PSFe1,PSFe2,TPSF,SNR)
	t2 = time.time()
	print "It took %.0f seconds to read in from the pickled files" %(t2-t1)	



if str(ZBlo) == "None":
    ZBlabel = 'ZBcutNone'
elif float(ZBlo)<float(ZBhi):
	print "Making the ZBcut in the range %s to %s" %(ZBlo, ZBhi)
	ZBlabel = 'ZBcut%s-%s' %(ZBlo,ZBhi)
	idx = np.where( np.logical_and(ZB>float(ZBlo), ZB<float(ZBhi)) )[0]
	ra=ra[idx]
	dec=dec[idx]
	e1=e1[idx]
	e2=e2[idx]
	w=w[idx]
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
# Save the additive bias - this isn't acually used for anythin atm.
np.savetxt('%s/Catalogues/K%s.Blind%s.ccorr_e1_e2.%s.dat'%(DIRECT,NorS,Blind,ZBlabel), np.c_[ np.mean(e1),np.mean(e2) ])


def Bootstrap_Error(nboot, samples, weights):
	N = len(samples)
	bt_samples = np.zeros(nboot)		 		# Will store mean of nboot resamples
	for i in range(nboot):
		idx = np.random.randint(0,N,N)			# Picks N random indicies with replacement
		bt_samples[i] = np.sum( weights[idx]*samples[idx] ) / np.sum( weights[idx] )
	return np.std( bt_samples )






def MeanQ_VS_XY(Q, X,Y):
	num_XY_bins = 20
	AvQ_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q,statistic='mean', bins=num_XY_bins)#, range=[[0.0,21000.],[0.0,21000.0]])
	#AvQ_grid[~np.isfinite(AvQ_grid)] = 0.
	return AvQ_grid

def Plot_XY_Grids(Q, lowlim, upplim, label, savename):
	cmap = plt.cm.rainbow
	cmap.set_bad('w', 1.)

	plt.figure(figsize=(8,7))
	plt.imshow(Q, cmap=cmap, interpolation='none', vmin=lowlim, vmax=upplim, origin='lower')
	plt.ylabel(r'$Y_{\rm{pos}}$')
	plt.xlabel(r'$X_{\rm{pos}}$')
	plt.colorbar(label=label)
	plt.savefig(savename)
	#plt.show()
	return



# Mean-e in bins VS Xpos,Ypos
e1_grid = MeanQ_VS_XY(e1, Xpos,Ypos)
e2_grid = MeanQ_VS_XY(e2, Xpos,Ypos)
#Plot_XY_Grids(1e5*e1_grid, -1, 1, r'$\langle e_1 \rangle\times 10^{-5}$', '%s/K%s.Blind%s.Meane1_VS_XY.%s.png' %(DIRECT,NorS,Blind,ZBlabel) )
#Plot_XY_Grids(1e5*e2_grid, -1, 1, r'$\langle e_2 \rangle\times 10^{-5}$', '%s/K%s.Blind%s.Meane2_VS_XY.%s.png' %(DIRECT,NorS,Blind,ZBlabel) )


# I also want to plot the PSF ellipticity and PSF residual on the chip, but the latter depends on the model fit
# so you have to read this in separately, from the file, made over on cuillin, used to make rho stats.
# The file you are reading in here contains the PSF measured from the stars, whereas the pickled file read in
# above is the PSF interpolated to the galaxy positions. This is why this data set is only 7% of the size of the former.
# So new variables needed I suppose. Add '_d2' to indicate it's from datafile2
PSFe1_d2_N, PSFe2_d2_N, delta_PSFe1_d2_N, delta_PSFe2_d2_N, Xpos_d2_N, Ypos_d2_N = np.load('PSF_rho_stats/LFver321/Catalogues/PSF_Data_N.npy').transpose()[[2,3,4,5,8,9],:]
PSFe1_d2_S, PSFe2_d2_S, delta_PSFe1_d2_S, delta_PSFe2_d2_S, Xpos_d2_S, Ypos_d2_S = np.load('PSF_rho_stats/LFver321/Catalogues/PSF_Data_S.npy').transpose()[[2,3,4,5,8,9],:]
# Append N and S.
PSFe1_d2 = np.append(PSFe1_d2_N, PSFe1_d2_S) 
PSFe2_d2 = np.append(PSFe2_d2_N, PSFe2_d2_S) 
delta_PSFe1_d2 = np.append(delta_PSFe1_d2_N, delta_PSFe1_d2_S)
delta_PSFe2_d2 = np.append(delta_PSFe2_d2_N, delta_PSFe2_d2_S)
Xpos_d2 = np.append(Xpos_d2_N, Xpos_d2_S)
Ypos_d2 = np.append(Ypos_d2_N, Ypos_d2_S)

PSFe1_grid = MeanQ_VS_XY(PSFe1_d2, Xpos_d2,Ypos_d2)
PSFe2_grid = MeanQ_VS_XY(PSFe2_d2, Xpos_d2,Ypos_d2)
delta_PSFe1_grid = MeanQ_VS_XY(delta_PSFe1_d2, Xpos_d2,Ypos_d2)
delta_PSFe2_grid = MeanQ_VS_XY(delta_PSFe2_d2, Xpos_d2,Ypos_d2)
Plot_XY_Grids(1e2*PSFe1_grid, -1, 1, r'$\langle e_{{\rm PSF},1} \rangle \times 10^{-2}$', 
				'%s/GeneralPlots/K%s.Blind%s.MeanPSFe1_VS_XY.png' %(DIRECT,NorS,Blind) )
Plot_XY_Grids(1e2*PSFe2_grid, -1, 1, r'$\langle e_{{\rm PSF},2} \rangle \times 10^{-2}$', 
				'%s/GeneralPlots/K%s.Blind%s.MeanPSFe2_VS_XY.png' %(DIRECT,NorS,Blind) )

Plot_XY_Grids(1e3*delta_PSFe1_grid, -1, 1, r'$\langle \delta e_{{\rm PSF},1} \rangle \times 10^{-3}$', 
				'%s/GeneralPlots/K%s.Blind%s.MeandeltaPSFe1_VS_XY.png' %(DIRECT,NorS,Blind) )
Plot_XY_Grids(1e3*delta_PSFe2_grid, -1, 1, r'$\langle \delta e_{{\rm PSF},2} \rangle \times 10^{-3}$', 
				'%s/GeneralPlots/K%s.Blind%s.MeandeltaPSFe2_VS_XY.png' %(DIRECT,NorS,Blind) )



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

def Plot_BinQx_VS_BinQy(Qx, Qy, weights, num_bins, labels, xlabel, ylabel, xlimits, ylimits, savename, Bootstrap):

	colors = ['magenta', 'dimgrey', 'darkblue', 'orange', 'lawngreen', 'red']
	f, ((ax1)) = plt.subplots(1, 1, figsize=(8,7))

	if len(Qx.shape) == 1:
		Qx = np.reshape(Qx, (1,len(Qx)))
		Qy = np.reshape(Qy, (1,len(Qy)))
    
	offset = [-0.01,0.01, -0.03,0.03]   

	# Average the corresponding Y's to the X's
	binx_centres = np.zeros([Qy.shape[0], num_bins])
	biny_centres = np.zeros([Qy.shape[0], num_bins])
	biny_centres_err = np.zeros([Qy.shape[0], num_bins])
	nboot = 10	# The number of bootstrap resamples to take
						# in estimating the errors.
	probability = np.linspace(0,1,num_bins+1)				
	for i in range(Qy.shape[0]):			
		binx_edges = mquantiles(Qx[i,:], prob=probability)				# Edges of the num_bins bins, each containing same-% of data points.
		
		for j in range(num_bins):
			idx = np.where(np.logical_and(Qx[i,:]>=binx_edges[j], Qx[i,:]<binx_edges[j+1]))[0]
			# print len(idx) / float(len(Qx[i,:]))
			binx_centres[i,j] = np.sum( weights[idx]*Qx[i,idx] ) / np.sum( weights[idx] ) 
			biny_centres[i,j] = np.sum( weights[idx]*Qy[i,idx] ) / np.sum( weights[idx] ) 
			if Bootstrap:
				biny_centres_err[i,j] = Bootstrap_Error(nboot, Qy[i,idx], weights[idx])
   			else:
				biny_centres_err[i,j] = 0.
 
		plt.errorbar(binx_centres[i,:], biny_centres[i,:], yerr=biny_centres_err[i,:], fmt='o', color=colors[i],label=labels[i])
		# +offset[i]*np.mean(binx_centres)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if ylimits != None:
		plt.ylim(ylimits)
	if xlimits != None:
		plt.plot( xlimits, [0.,0.], 'k--' )
		plt.xlim(xlimits)
	plt.legend(loc='best', frameon=False)
	plt.savefig(savename)
	#plt.show()
	return np.dstack(( binx_centres, biny_centres, biny_centres_err ))



# ------------------------------------------------------------- 1-POINT STATISTICS PLOTS ------------------------------------------------------------- #
# NOTE ON RETURNED ARRAY: bin_x_y_yerr[:,:,i] --> i=0 for x, i=1 for y, i=2 for yerr.

print "Exiting as the next 1pt plots are expensive (~100s) to produce due to bootstrapping errors."
sys.exit()


Bootstrap = True	# Bootstrap errors? Takes ~nboot times longer! (nboot defined in Plot_BinQx_VS_BinQy)

# MEAN ELLIPTICITY VS X/Ypos
#bin_Xpos_e_eerr  = Plot_BinQx_VS_BinQy(np.vstack((Xpos,Xpos)), np.vstack((e1, e2)), w, 10, 
#[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$X_{\rm{pos}}$', r'$\langle e \rangle$', 
#None, '%s/GeneralPlots/K%s.Blind%s.Xpos_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )

#bin_Ypos_e_eerr = Plot_BinQx_VS_BinQy(np.vstack((Ypos,Ypos)), np.vstack((e1, e2)), w, 10, 
#[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$Y_{\rm{pos}}$', r'$\langle e \rangle$', 
#None, '%s/GeneralPlots/K%s.Blind%s.Ypos_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )

t1 = time.time()
# MEAN ELLIPTICITY VS PSF ELLIPTICITY
bin_PSFe_e_eerr = Plot_BinQx_VS_BinQy(np.vstack((PSFe1,PSFe2)), 1000*np.vstack((e1, e2)), w, 20, 
[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$e_{\rm{PSF}}$', r'$\langle e \rangle \times 10^{-3}$',
[-0.04,0.04], [-2.5,2.5], 
'%s/GeneralPlots/K%s.Blind%s.ePSF_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap )
t2 = time.time()
print "It took %.0f s to bin and plot e VS PSFe" %(t2-t1)


# MEAN ELLIPTICITY VS PSF SIZE
bin_TPSF_e_err = Plot_BinQx_VS_BinQy(np.vstack((TPSF,TPSF)), 1000*np.vstack((e1, e2)), w, 20, 
[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$T_{\rm{PSF}}$', r'$\langle e \rangle \times 10^{-3}$',
[2.5,5.5], [-2.5,2.5], 
'%s/GeneralPlots/K%s.Blind%s.TPSF_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )
print "Finished binning and plotting e VS TPSF" 

# MEAN ELLIPTICITY VS MAGNITUDE
bin_MAG_e_err = Plot_BinQx_VS_BinQy(np.vstack((MAG,MAG)), 1000*np.vstack((e1, e2)), w, 20, 
[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'Magnitude', r'$\langle e \rangle \times 10^{-3}$',
[20,25.5], [-2.5,2.5], 
'%s/GeneralPlots/K%s.Blind%s.MAG_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), Bootstrap  )
print "Finished binning and plotting e VS MAG" 

# MEAN ELLIPTICITY VS ZB
bin_ZB_e_err = Plot_BinQx_VS_BinQy(np.vstack((ZB,ZB)), 1000*np.vstack((e1, e2)), w, 5, 
[r'$\langle e_1 \rangle$', r'$\langle e_2 \rangle$'], r'$z_{\rm B}$', r'$\langle e \rangle \times 10^{-3}$',
[0.1,1.5], [-2.5,2.5], 
'%s/GeneralPlots/K%s.Blind%s.ZB_VS_e.%s.png' %(DIRECT,NorS,Blind,ZBlabel), False  )
print "Finished binning and plotting e VS ZB" 
plt.show()























