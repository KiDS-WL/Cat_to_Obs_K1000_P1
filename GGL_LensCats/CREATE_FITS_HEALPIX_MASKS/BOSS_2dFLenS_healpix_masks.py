

    ##############################
    ##  MFP_K1000.py            ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.03.25      ##
    ##############################


import os
import os.path as osp
import time
import subprocess as spc

import numpy as np
import scipy as sp
import astropy.io.fits as fits
import healpy as hp
import treecorr as tree

import commonFunctions as cf
import HEALPixFunctions as hpf


################################################################################
## Parameters

class Parameters:
  KiDSPath    = 'data/KiDS/'
  dataPath    = 'data/mockFootprint/'
  absDataPath = '/disk05/calin/91_Data/mockFootprint/'
  
  ## Mask parameters
  area_BOSS          = 9329 ## [deg^2]
  area_BOSS_reduced  = 1274.319868 ## From my own calculations
  area_BOSS_wcs      = 408.321
  area_BOSS_4Band    = 339.298
  area_BOSS_9Band    = 319.506
  
  area_2dFLenS_SGP   = 510.803964 ## [deg^2]
  area_2dFLenS_wcs   = 424.508017
  area_2dFLenS_gri   = 355.283139
  area_2dFLenS_9Band = 341.888289
  
  area_KiDS          = 773.286 ## [deg^2]
  area_KiDS_North    = 334.138
  area_KiDS_South    = 439.148
  
  area_KiDS_North_new = 371.801
  area_KiDS_South_new = 401.485
  
  ## Galaxy number density
  n_gal_BOSS_reduced_z0 = 0.014496
  n_gal_BOSS_reduced_z1 = 0.016595
  
  n_gal_BOSS_wcs_z0     = 0.014437
  n_gal_BOSS_wcs_z1     = 0.016265
  
  n_gal_2dFLenS_SGP_z0  = 0.005813
  n_gal_2dFLenS_SGP_z1  = 0.006067
  
  n_gal_2dFLenS_wcs_z0  = 0.005857
  n_gal_2dFLenS_wcs_z1  = 0.006031
  
  n_gal_2dFLenS_gri_z0  = 0.002891
  n_gal_2dFLenS_gri_z1  = 0.003677
  
################################################################################
## Functions related to masks - I

## This function load BOSS random catalogues
def loadFitsLenCat(surveyTag, zInd, bitMaskTag='reduced'):
  P = Parameters()
  if bitMaskTag in ['all', 'reduced', 'SGP']: ## No selection
    bitMask = 000000
  elif bitMaskTag == 'wcs': ## KiDS wcs
    bitMask = 0x4000
  elif bitMaskTag == 'gri':
    bitMask = 0x6FFC ## KiDS gri overlap
  elif bitMaskTag == '9Band':
    bitMask = 0x681C ## KiDS 9-band overlap
  else:
    raise ValueError('Bad bit mask option: \"%s\"' % bitMaskTag)
  
  name = '%sKiDS-1000_GGLCATS/%s_z%d.fits' % (P.KiDSPath, surveyTag, zInd+1)
  data = fits.getdata(name, 1)
  print('Loaded \"%s\"' % name)
  
  flag = data.field('KIDSMASK')
  ind  = np.logical_not(np.array(flag.astype(int) & bitMask, dtype=bool))
  return data[ind]

## This function loads BOSS random catalogues & pour them onto a HEALPix map. 
def saveFitsCountMap_BOSS(nside, bitMaskTag='wcs'):
  P = Parameters()
  nbPix = 12 * nside * nside
  full  = np.zeros(nbPix, dtype=int)
  
  ## Fill catalogues
  for zInd in range(2):
    data = loadFitsLenCat('BOSS_random', zInd, bitMaskTag=bitMaskTag)
    RA   = data.field('ALPHA_J2000')
    DEC  = data.field('DELTA_J2000')
    pix  = hpf.RADECToPatch(nside, RA, DEC)
    for i in pix:
      full[i] += 1
  
  ## Save
  name = '%sKiDS-1000_for_mocks/countMap_BOSS_%s_nside%d.fits' % (P.KiDSPath, bitMaskTag, nside)
  hpf.saveFitsFullMap(name, full, verbose=True)
  return

def saveFitsCountMap_overlap(surveyTag_K, surveyTag_L, nside_L):
  P = Parameters()
  nside_K = 4096
  
  name = '%sKiDS-1000_for_mocks/countMap_%s_nside%d.fits' % (P.KiDSPath, surveyTag_L, nside_L)
  count_L = hpf.loadFitsFullMap(name)
  count_L = hpf.increaseResolution(count_L, nside_K)
  
  name = '%sKiDS-1000_for_mocks/mask_%s_fromArea_nside%d.fits' % (P.KiDSPath, surveyTag_K, nside_K)
  mask_K = hpf.loadFitsFullMap(name)
  
  ind = mask_K.astype(bool)
  del mask_K
  
  count_L[~ind] = 0
  del ind
  
  ## Save
  surveyTag_o = 'BOSS_KiDS_overlap' if 'BOSS' in surveyTag_L else '2dFLenS_KiDS_overlap'
  name = '%sKiDS-1000_for_mocks/countMap_%s_nside%d.fits' % (P.KiDSPath, surveyTag_o, nside_K)
  hpf.saveFitsFullMap(name, count_L)
  del count_L
  return

## 'BOSS_wcs' is called
def saveFitsMask_fromCountMap(surveyTag):
  P = Parameters()
  if surveyTag == 'BOSS_reduced':
    nside = 2048
  elif surveyTag == 'BOSS_wcs':
    nside = 2048
  elif surveyTag == '2dFLenS_SGP':
    nside = 4096
  elif surveyTag == '2dFLenS_wcs':
    nside = 4096
  else:
    raise NotImplementedError('surveyTag = \"%s\" not implemented' % surveyTag)
  
  name = '%sKiDS-1000_for_mocks/countMap_%s_nside%d.fits' % (P.KiDSPath, surveyTag, nside)
  mask = hpf.loadFitsFullMap(name)
  mask = np.fmin(mask, 1)
  
  if nside == 2048:
    nside2 = 4096
    mask = hpf.increaseResolution(mask, nside2)
    name = '%sKiDS-1000_for_mocks/mask_%s_fromCountMap2048_nside%d.fits' % (P.KiDSPath, surveyTag, nside2)
    hpf.saveFitsFullMap(name, mask)
    return
  
  ## Save
  name = '%sKiDS-1000_for_mocks/mask_%s_fromCountMap_nside%d.fits' % (P.KiDSPath, surveyTag, nside)
  hpf.saveFitsFullMap(name, mask)
  return

# This function combines the 2dFLenS mask and BOSS mask into one
def saveFitsLensMask():
    P = Parameters()
    name   = '%sKiDS-1000_for_mocks/mask_BOSS_wcs_fromCountMap2048_nside4096.fits' % P.KiDSPath
    mask_B = hpf.loadFitsFullMap(name)
    name   = '%sKiDS-1000_for_mocks/mask_2dFLenS_wcs_fromCountMap_nside4096.fits' % P.KiDSPath
    mask_2 = hpf.loadFitsFullMap(name)
                        
    mask_L = mask_B + mask_2
    mask_L = np.fmin(mask_L, 1)
                                
    name   = '%sKiDS-1000_for_mocks/mask_BOSS_2dFLenS_wcs_nside4096.fits' % P.KiDSPath
    hpf.saveFitsFullMap(name, mask_L)
    return

## Then I called the following & used the output of the 2nd line
##  saveFitsCountMap_BOSS(2048, 'wcs') ## Need external
##  saveFitsMask_fromCountMap('BOSS_wcs')

###############################################################################

