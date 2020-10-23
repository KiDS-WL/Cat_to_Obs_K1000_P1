#NZDATA  = T          // This sentinel marks the extension as n(z) data
#EXTNAME = NAME       // The name of this n(z) kernel.
#NBIN    = 5          // Integer number of tomographic bins
#NZ      = 100        // Integer number of histogram bins
#The extension must then contain these data columns:

#Z_LOW   8-byte real  // Real, the z value for the lower end of each redshift histogram bin
#Z_MID   8-byte real  // Real, the z value for a middle point of each redshift histogram bin
#Z_HIGH  8-byte real  // Real, the z value for the upper end of each redshift histogram bin
#BIN1    8-byte real  // Real, the n(z) value for this histogram bin for the first tomographic bin
#etc.    BIN2, BIN3, etc.


import sys
import os
import numpy as np
import pylab as plt
from   matplotlib.ticker import ScalarFormatter
import pyfits


# KiDS-1000 values for blindA
nBins_KIDS=5
nBins_BOSS=2
neff=[0.5696062756935643, 1.1420051774183186, 1.820071990524645, 1.257355966105446, 1.3150376529887304]



def MakeNofz_fits(input_files,outputfileName,OneBin_nofzFileName,neff=[],single_bin=True,type='lowerEdge',suffix='SAMPLE'):
	nBins=len(input_files)
	print('I got '+str(nBins)+' files from input. Type is set to '+type)
	cols=[]
	for bin1 in range(nBins):
		file= open(input_files[bin1])
		nofZ=np.loadtxt(file,comments='#')
		if(bin1==0):
			if(single_bin):
				z_vec=nofZ[:,1]*neff[bin1]
			DeltaZ=nofZ[1,0]-nofZ[0,0]
			if(type=='lowerEdge'):
				Z_LOW=nofZ[:,0]
				Z_HIGH=nofZ[:,0]+DeltaZ
				Z_MID=Z_LOW+DeltaZ/2.
			elif(type=='middle'):
				Z_MID=nofZ[:,0]
				Z_LOW=nofZ[:,0]-DeltaZ/2.
				Z_HIGH=nofZ[:,0]+DeltaZ/2.
			elif(type=='upperEdge'):
				Z_HIGH=nofZ[:,0]
				Z_MID=nofZ[:,0]-DeltaZ/2.
				Z_LOW=nofZ[:,0]-DeltaZ
			else:
				print('not a recognised bin type, exiting now ...')
				exit(1)

			cols.append(pyfits.Column(name='Z_lOW', format='D', array=Z_LOW))
			cols.append(pyfits.Column(name='Z_HIGH', format='D', array=Z_HIGH))
			cols.append(pyfits.Column(name='Z_MID', format='D', array=Z_MID))
		else:
			if(single_bin):
				z_vec+=nofZ[:,1]*neff[bin1]
		cols.append(pyfits.Column(name='BIN'+str(bin1+1), format='D', array=nofZ[:,1]))

	new_cols = pyfits.ColDefs(cols)
	#what happened here?
	#if python version older than 3.3
	#hdulist_new = pyfits.new_table(data.columns+new_cols)
	#else
	hdulist_new = pyfits.BinTableHDU.from_columns(new_cols)
	hdulist_new.header['NZDATA'] = True
	hdulist_new.header['EXTNAME'] = 'NZ_'+suffix
	hdulist_new.header['NBIN'] = 5
	hdulist_new.header['NZ'] = len(Z_LOW)
	hdulist_new.writeto(outputfileName)
	# now one bin
	if(single_bin):
		cols = [] 
		cols.append(pyfits.Column(name='Z_lOW', format='D', array=Z_LOW))
		cols.append(pyfits.Column(name='Z_HIGH', format='D', array=Z_HIGH))
		cols.append(pyfits.Column(name='Z_MID', format='D', array=Z_MID))
		cols.append(pyfits.Column(name='BIN1', format='D', array=z_vec))
		new_cols = pyfits.ColDefs(cols)
		#what happened here?
		#if python version older than 3.3
		#hdulist_new = pyfits.new_table(data.columns+new_cols)
		#else
		#OneBin_nofzFileName='Version2/Nz_DIR/Nz_DIR_Mean/nofZ1bin.fits'
		hdulist_new = pyfits.BinTableHDU.from_columns(new_cols)
		hdulist_new.header['NZDATA'] = True
		hdulist_new.header['EXTNAME'] = 'NZ_'+suffix
		hdulist_new.header['NBIN'] = 1
		hdulist_new.header['NZ'] = len(Z_LOW)
		outputfileName=OneBin_nofzFileName
		hdulist_new.writeto(outputfileName)

cat_version_out = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'
cat_version_in  = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_Fid'
for blind in 'A','B','C':
	name_in='K1000_NS_'+cat_version_in+'_blind'+blind
	name_out = 'K1000_NS_'+cat_version_out+'_blind'+blind
	OutputFileName='SOM_N_of_Z/'+name_out+'.fits'
	OneBin_nofzFileName='SOM_N_of_Z/'+name_out+'_single_bin.fits'
	input_files=[]
	for bin1 in range(nBins_KIDS):
		fileNameInput='SOM_N_of_Z/'+name_in+'_TOMO'+str(bin1+1)+'_Nz.asc'
		input_files.append(fileNameInput)
# 
	MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,single_bin=False,type='lowerEdge',suffix='source')

# BOSS inputs
nBins_BOSS=2
input_files=[]
name_in="BOSS_and_2dFLenS_n_of_z"
OutputFileName = '../../boss/nofz/'+name_in + '.fits'
OneBin_nofzFileName = '../../boss/nofz/'+name_in+'_single_bin.fits'
res=0.01
z_arr=np.arange(0,0.8,res)
for bin1 in range(nBins_BOSS):
	fileNameInput='../../boss/nofz/'+name_in+str(bin1+1)+'_res_'+ '%0.2f' % res+'.txt'
	file= open(fileNameInput)
	nofZ_in=np.loadtxt(file,comments='#')
	nofZ_out=np.zeros((len(z_arr),2))
	nz_start = nofZ_in[0,0]
	nz_end   = nofZ_in[-1,0]
	nofZ_out[:,0]=z_arr
	start_arg=np.argwhere(z_arr==nz_start)[0,0]
	end_arg  =np.argwhere(z_arr==nz_end)[0,0]
	nofZ_out[start_arg:end_arg+1,1]=nofZ_in[:,1]
	fileNameOut='../../boss/nofz/'+name_in+str(bin1+1)+'_res_'+ '%0.2f' % res+'_extended.txt'
	np.savetxt(fileNameOut,nofZ_out,fmt='%0.3f')
	input_files.append(fileNameOut)

neff=[0.02,0.02]
MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,single_bin=False,type='lowerEdge',suffix='lens')


