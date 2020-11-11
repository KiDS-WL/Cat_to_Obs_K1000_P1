# 30/03/2018, B. Giblin, Postdoc., Edinburgh,
# Edited from original code by M. Asgari.
# Read in either BOSS (lens) or K1000 (Source) N(z)
# and make a single FITS file used by kcap ini file, magnifcation_gt.ini 

import sys
import os
import numpy as np
import pylab as plt
from   matplotlib.ticker import ScalarFormatter
import pyfits

DATA_TYPE = "SOURCE"                   # If "SOURCE", it's K1000
                                       # If "LENS", it's determined by LENS_TYPE
OL_Tag = '_OutlierPeaksInBins12345'    # If there's high-z outliers in the source nofz. ('' if not).
LENS_TYPE = "GAMA"

                     
# Use this z-axis for everything                     
z = np.linspace( 0., 2.0, 201 )
                     
def MakeNofz_fits(input_files,outputfileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix='SAMPLE'):
        nBins=len(input_files)
        print('I got '+str(nBins)+' files from input. Type is set to '+type)
        cols=[]
        for bin1 in range(nBins):
                file= open(input_files[bin1])
                z_tmp,dens_tmp=np.loadtxt(file, usecols=(0,1), unpack=True)
                density = np.interp( z, z_tmp, dens_tmp )
                
                # We want the lenses to have the same z-axis as the K1000 sources
                if DATA_TYPE == "LENS":
                        print("Saving a LENS nofz - copying K1000 z-axis array.")
                        K1000_file = 'examples/SOM_NofZ/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_Fid_blindA.fits'
                        f = pyfits.open(K1000_file)
                        z = f[1].data['Z_LOW']
                        density = np.interp( z, z_tmp, dens_tmp)
                else:
                        z = np.copy( z_tmp )
                        density = np.copy( dens_tmp )
                        
                
                if(bin1==0):
                        density_all=density*neff[bin1]
                        DeltaZ=z[1]-z[0]
                        if(type=='lowerEdge'):
                                Z_LOW=z
                                Z_HIGH=z+DeltaZ
                                Z_MID=z+DeltaZ/2.
                        elif(type=='middle'):
                                Z_MID=z
                                Z_LOW=z-DeltaZ/2.
                                Z_HIGH=z+DeltaZ/2.
                        elif(type=='upperEdge'):
                                Z_HIGH=z
                                Z_MID=z-DeltaZ/2.
                                Z_LOW=z-DeltaZ
                        else:
                                print('not a recognised bin type, exiting now ...')
                                exit(1)

                        cols.append(pyfits.Column(name='Z_lOW', format='D', array=Z_LOW))
                        cols.append(pyfits.Column(name='Z_HIGH', format='D', array=Z_HIGH))
                        cols.append(pyfits.Column(name='Z_MID', format='D', array=Z_MID))
                else:
                        density_all+=density*neff[bin1]
                cols.append(pyfits.Column(name='BIN'+str(bin1+1), format='D', array=density))

        new_cols = pyfits.ColDefs(cols)
        hdulist_new = pyfits.BinTableHDU.from_columns(new_cols)
        hdulist_new.header['NZDATA'] = True
        hdulist_new.header['EXTNAME'] = 'NZ_'+suffix
        hdulist_new.header['NBIN'] = 5
        hdulist_new.header['NZ'] = len(Z_LOW)
        hdulist_new.writeto(outputfileName, clobber=True)
	# now one bin
        cols = [] 
        cols.append(pyfits.Column(name='Z_lOW', format='D', array=Z_LOW))
        cols.append(pyfits.Column(name='Z_HIGH', format='D', array=Z_HIGH))
        cols.append(pyfits.Column(name='Z_MID', format='D', array=Z_MID))
        cols.append(pyfits.Column(name='BIN1', format='D', array=density_all))
        new_cols = pyfits.ColDefs(cols)
        hdulist_new = pyfits.BinTableHDU.from_columns(new_cols)
        hdulist_new.header['NZDATA'] = True
        hdulist_new.header['EXTNAME'] = 'NZ_'+suffix
        hdulist_new.header['NBIN'] = 1
        hdulist_new.header['NZ'] = len(Z_LOW)
        outputfileName=OneBin_nofzFileName
        hdulist_new.writeto(outputfileName, clobber=True)

if DATA_TYPE == "SOURCE":
        nBins=5
        #neff=[0.5696062756935643, 1.1420051774183186, 1.820071990524645, 1.257355966105446, 1.3150376529887304]
        neff=[1.,1.,1.,1.,1.]
        for blind in 'A','B','C':
                inDIR = '/home/bengib/KiDS1000_NullTests/Cat_to_Obs_K1000_P1/Shear_Ratio_Test/SOURCECATS/SOM_NofZ'
                filename = 'K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_Fid_blind'+blind

                outDIR = 'examples/SOM_NofZ/'
                OutputFileName = outDIR + filename + OL_Tag + '.fits'
                OneBin_nofzFileName = outDIR + filename + OL_Tag + '_single_bin.fits'
                input_files=[]
                for bin1 in range(nBins):
                        fileNameInput = inDIR + filename +'_TOMO'+str(bin1+1)+'_Nz' +OL_Tag+'.asc'
                        input_files.append(fileNameInput)
                MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix=DATA_TYPE)

elif DATA_TYPE == "LENS":
        nBins=5
        neff=[1., 1., 1., 1., 1.]
        inDIR = '/home/bengib/KiDS1000_NullTests/Cat_to_Obs_K1000_P1/Shear_Ratio_Test/LENSCATS/%s_data/' %LENS_TYPE
        outDIR = 'examples/%s_Nofz/' %LENS_TYPE
        filename = 'nofz_cat_5Z'
        OutputFileName = outDIR + filename + '.fits'
        OneBin_nofzFileName = outDIR + filename + '_single_bin.fits' 

        input_files=[]
        for bin1 in range(nBins):
                fileNameInput = inDIR + filename + '_%s.dat' %(bin1+1)
                input_files.append(fileNameInput)
        MakeNofz_fits(input_files,OutputFileName,OneBin_nofzFileName,neff,type='lowerEdge',suffix=DATA_TYPE)
