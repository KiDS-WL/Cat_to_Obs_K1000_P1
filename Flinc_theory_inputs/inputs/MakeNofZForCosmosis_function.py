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


def MakeNofz_fits(input_files,outputfileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix='SAMPLE'):
	nBins=len(input_files)
	print('I got '+str(nBins)+' files from input. Type is set to '+type)
	cols=[]
	for bin1 in range(nBins):
		file= open(input_files[bin1])
		nofZ=np.loadtxt(file,comments='#')
		if(bin1==0):
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


nBins_lens=2
# neff=[1.5099137858524687,1.5413199175708638 , 1.6196461014383043, 0.8268988140297754]
neff=[1,1]
OutputFileName='lens_flinc_nofz.fits'
OneBin_nofzFileName='lens_flinc_nofz_1bin.fits'
input_files=[]
for bin1 in range(nBins_lens):
	fileNameInput='numida/lens/'+'nOfZ_hist_BOSSA_tomo'+str(bin1)+'.dat'
	input_files.append(fileNameInput)


MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix='lens_flinc')


nBins_source=5
# neff=[1.5099137858524687,1.5413199175708638 , 1.6196461014383043, 0.8268988140297754]
neff=[1,1,1,1,1]
OutputFileName='source_flinc_nofz.fits'
OneBin_nofzFileName='source_flinc_nofz_1bin.fits'
input_files=[]
for bin1 in range(nBins_source):
	fileNameInput='numida/source/'+'nOfZ_hist_KiDSVDA_tomo'+str(bin1)+'.dat'
	input_files.append(fileNameInput)


MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix='source_flinc')
