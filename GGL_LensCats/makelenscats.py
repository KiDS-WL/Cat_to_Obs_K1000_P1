########################################################################
# Code to generate BOSS and 2dFLenS data and random lens catalogues    #
# in the KiDS regions, including magnitude weights, from the           #
# publicly-available datasets. 
# Author:  Chris Blake
# Questions to: cblake@swin.edu.au
# Original version 13th May 2019                                        #
# History
#
#
########################################################################

import sys
import numpy as np
import scipy.spatial
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

#============================================
# Generate the lens catalogues
def makecats(ired):

# Read in KiDS photometric catalogue both N and S
# (rmin_N,rmax_N,dmin_N,dmax_N) are returned as the precise boundaries of the KiDS-N region
# (rmin_S,rmax_S,dmin_S,dmax_S) are returned as the precise boundaries of the KiDS-S region
  raskids,deckids,grcolkids,ricolkids,rmagkids,rmin_N,rmax_N,dmin_N,dmax_N,rmin_S,rmax_S,dmin_S,dmax_S = readkids()

# Read in BOSS data and random lenses

# K1000-N has rough min/max RA (238.6, 128.4) and min/max DEC (-4.1, 3.1)
# We want to measure GGL out to 6 degrees so we should select the BOSS region
# to be 6.0 degrees larger than the KiDS region (buffer for good measure)

# (R.A., Dec.) boundaries to use for BOSS catalogue
  edge = 6.0
  rmin,rmax,dmin,dmax = rmin_N-edge,rmax_N+edge,dmin_N-edge,dmax_N+edge
  rasbossdat,decbossdat,redbossdat,weicompbossdat,weifkpbossdat,nbossdat = readboss(1,ired,rmin,rmax,dmin,dmax)
  rasbossran,decbossran,redbossran,weicompbossran,weifkpbossran,nbossran = readboss(2,ired,rmin,rmax,dmin,dmax)
  
# You might want to sub-sample the BOSS randoms to increase your speed
# We need high signal-to-noise random gamma_t signals though so we do not sub-sample here
# Sub-sample BOSS randoms to 40x data for consistency with 2dFLenS
#  cut = np.random.choice(nbossran,40*nbossdat,replace=False)
#  rasbossran,decbossran,redbossran,weicompbossran,weifkpbossran = rasbossran[cut],decbossran[cut], \
#                                                                  redbossran[cut],weicompbossran[cut],weifkpbossran[cut]
#  nbossran = 40*nbossdat
#  print 'Sub-sampled BOSS randoms to',nbossran,'lenses'

# Read in 2dFLenS data and random lenses 
# We do not need to apply ra/dec cuts here as 2dFLenS is designed to overlap with KiDS
  ras2dfdat,dec2dfdat,red2dfdat,weifkp2dfdat,n2dfdat = read2dflens(1,ired)
  ras2dfran,dec2dfran,red2dfran,weifkp2dfran,n2dfran = read2dflens(2,ired)
# completeness weights=1 for 2dFLenS
  weicomp2dfdat,weicomp2dfran = np.ones(n2dfdat),np.ones(n2dfran)

# Find KiDS mags/colours of BOSS and 2dFLenS data lenses
# iphot is a flag which is 1 if there is a photometry match, otherwise 0

  grcolbossdat,ricolbossdat,rmagbossdat,iphotbossdat = matchlenstosource(rasbossdat,decbossdat,raskids,deckids,\
                                                                         grcolkids,ricolkids,rmagkids,\
                                                                         rmin_N,rmax_N,dmin_N,dmax_N,\
                                                                         rmin_S,rmax_S,dmin_S,dmax_S)
  grcol2dfdat,ricol2dfdat,rmag2dfdat,iphot2dfdat = matchlenstosource(ras2dfdat,dec2dfdat,raskids,deckids,\
                                                                     grcolkids,ricolkids,rmagkids,\
                                                                     rmin_N,rmax_N,dmin_N,dmax_N,\
                                                                     rmin_S,rmax_S,dmin_S,dmax_S)

# Determine magnitude weights of 2dFLenS data with BOSS as reference
  cutboss = (iphotbossdat > 0)
  cut2df = (iphot2dfdat > 0)
  magsbossdat = np.dstack([grcolbossdat[cutboss],ricolbossdat[cutboss],rmagbossdat[cutboss]])[0]
  mags2dfdat = np.dstack([grcol2dfdat[cut2df],ricol2dfdat[cut2df],rmag2dfdat[cut2df]])[0]
  weimag2dfdat = np.ones(n2dfdat)
  weimag2dfdat[cut2df] = calcmagweights(mags2dfdat,magsbossdat,weicompbossdat[cutboss])

# Determine magnitude weights of BOSS data with 2dFLenS as reference
  weimagbossdat = np.ones(nbossdat)
  weimagbossdat[cutboss] = calcmagweights(magsbossdat,mags2dfdat,weicomp2dfdat[cut2df])
# magnitudes=0 and weights=1 for randoms
  
  iphotbossran,weimagbossran,grcolbossran,ricolbossran,rmagbossran = np.zeros(nbossran,dtype='int'),np.ones(nbossran),\
                                                                     np.zeros(nbossran),np.zeros(nbossran),np.zeros(nbossran)
  iphot2dfran,weimag2dfran,grcol2dfran,ricol2dfran,rmag2dfran = np.zeros(n2dfran,dtype='int'),np.ones(n2dfran),\
                                                                np.zeros(n2dfran),np.zeros(n2dfran),np.zeros(n2dfran)
# Write out fits file catalogues
  outfile = OUTDIR +'/BOSS_data_z' + str(ired) + '.fits'
  writelenscat(outfile,rasbossdat,decbossdat,redbossdat,weicompbossdat,weifkpbossdat,\
                       iphotbossdat,weimagbossdat,grcolbossdat,ricolbossdat,rmagbossdat)
  outfile = OUTDIR +'/BOSS_random_z' + str(ired) + '.fits'
  writelenscat(outfile,rasbossran,decbossran,redbossran,weicompbossran,weifkpbossran,\
                       iphotbossran,weimagbossran,grcolbossran,ricolbossran,rmagbossran)
  outfile = OUTDIR +'/2dFLenS_data_z' + str(ired) + '.fits'
  writelenscat(outfile,ras2dfdat,dec2dfdat,red2dfdat,weicomp2dfdat,weifkp2dfdat,\
                       iphot2dfdat,weimag2dfdat,grcol2dfdat,ricol2dfdat,rmag2dfdat)
  outfile = OUTDIR +'/2dFLenS_random_z' + str(ired) + '.fits'
  writelenscat(outfile,ras2dfran,dec2dfran,red2dfran,weicomp2dfran,weifkp2dfran,\
                       iphot2dfran,weimag2dfran,grcol2dfran,ricol2dfran,rmag2dfran)
  return

#============================================
# Read in KiDS photometric catalogue
def readkids():
  print '\nReading in KiDS bright source data...'

  raskids,deckids,grcolkids,ricolkids,rmagkids = [],[],[],[],[]

  # Read in KiDS N and KiDS S and combine into a single data vector
  # Do not do this with anything other than a bright sample

  for ireg in range(1,3):
    if (ireg == 1):
      datfile = 'K1000_N_'+KiDS_VER
    else:
      datfile = 'K1000_S_'+KiDS_VER

    hdulist = fits.open(KiDS_DIR+'/'+datfile)

    # The KiDS catalogues are ldac tables in the 2nd extension
    table = hdulist[2].data
    
    raskids1 = table.field('ALPHA_J2000')
    # trick to deal with the zero-crossing
    if (ireg == 2):
      raskids1[raskids1 > 180.] -= 360.

    deckids1 = table.field('DELTA_J2000')
    raskids = np.append(raskids,raskids1)
    deckids = np.append(deckids,deckids1)
    gmaggaap = table.field('MAG_GAAP_g')
    rmaggaap = table.field('MAG_GAAP_r')
    imaggaap = table.field('MAG_GAAP_i')
    rmagtot = table.field('MAG_AUTO')
    grcolkids = np.append(grcolkids,gmaggaap-rmaggaap)
    ricolkids = np.append(ricolkids,rmaggaap-imaggaap)
    rmagkids = np.append(rmagkids,rmagtot)

    if (ireg == 1):
      rmin_N,rmax_N,dmin_N,dmax_N = np.amin(raskids1),np.amax(raskids1),np.amin(deckids1),np.amax(deckids1)
    else:
      rmin_S,rmax_S,dmin_S,dmax_S = np.amin(raskids1),np.amax(raskids1),np.amin(deckids1),np.amax(deckids1)
    hdulist.close()

  print len(raskids),'KiDS sources'
  # trick to deal with the zero-crossing
  raskids[raskids < 0.] += 360.
  
  return raskids,deckids,grcolkids,ricolkids,rmagkids,rmin_N,rmax_N,dmin_N,dmax_N,rmin_S,rmax_S,dmin_S,dmax_S

#============================================
# Read in BOSS lenses: datopt -- 1) data 2) random
def readboss(datopt,ired,rmin,rmax,dmin,dmax):

  # this will break if the file name is longer than 500
  datfile=np.chararray(2, itemsize=500)

  if (datopt == 1):
    print '\nReading in BOSS data lenses...'
  else:
    print '\nReading in BOSS random lenses...'
  if (ired == 1):
    zmin,zmax = 0.2,0.5
  elif (ired == 2):
    zmin,zmax = 0.5,0.75
  else:
    zmin,zmax = 0.4,0.6  # Overlap bin - not used
  if (datopt == 1):
    datfile[0] = 'galaxy_DR12v5_CMASSLOWZTOT_North.fits'
    nfiles = 1
  else:
    datfile[0] = 'random0_DR12v5_CMASSLOWZTOT_North.fits'
    datfile[1] = 'random1_DR12v5_CMASSLOWZTOT_North.fits'
    nfiles = 2

  #rasboss_out,decboss_out,redboss_out,weicompboss_out,weifkpboss_out = [],[],[],[],[]

  # read in files - there are two for the randoms
  for ifile in range(nfiles):
    hdulist = fits.open(BOSS_DIR+'/'+datfile[ifile])
    table = hdulist[1].data
    rasboss = table.field('RA')
    decboss = table.field('DEC')
    redboss = table.field('Z')
    weifkpboss = table.field('WEIGHT_FKP')
    if (datopt == 1):
      weicp = table.field('WEIGHT_CP')
      weinoz = table.field('WEIGHT_NOZ')
      weisys = table.field('WEIGHT_SYSTOT')
      weicompboss = weisys*(weinoz+weicp-1.)
      print len(rasboss),'BOSS lenses'
    else:
      weicompboss = np.ones(len(rasboss))
      print len(rasboss),'BOSS randoms', ifile
    hdulist.close()

    #ra/dec cuts
    cut = (rasboss > rmin) & (rasboss < rmax) & (decboss > dmin) & (decboss < dmax) & (redboss > zmin) & (redboss < zmax)

    rasboss_out=rasboss[cut]    
    decboss_out=decboss[cut]
    redboss_out=redboss[cut]
    weicompboss_out=weicompboss[cut]
    weifkpboss_out=weifkpboss[cut] 
  
    if (ifile>0):
      rasboss_out=np.append(rasboss_out,rasboss[cut]) 

  nboss = len(rasboss_out)
  print 'Cut to',nboss,'BOSS lenses with',rmin,'< R.A. <',rmax,dmin,'< Dec. <',dmax,zmin,'< z <',zmax
  return rasboss_out,decboss_out,redboss_out,weicompboss_out,weifkpboss_out,nboss

#============================================
# Read in 2dFLenS lenses: datopt -- 1) data 2) random
def read2dflens(datopt,ired):
  if (datopt == 1):
    print '\nReading in 2dFLenS data lenses...'
    nset = 1
  else:
    print '\nReading in 2dFLenS random lenses...'
    nset = 40
  ras2df,dec2df,red2df,weifkp2df = [],[],[],[]
  for iset in range(nset):
    for ireg in range(1,3):
      if (ireg == 1):
        creg = '_atlas_kidsn_160105'
      else:
        creg = '_atlas_kidss_160105'
      if (ired == 1):
        cred = '_bz1'
      elif (ired == 2):  # this is not a typo - Since the random catalogues have been created I have renamed 
                         # bins 2 and 3 as we will not use the overlap bin and it's easier to script over bins 1&2
                         # rather than over bins 1&3
        cred = '_bz3'
      else:
        cred = '_bz2'
      if (datopt == 1):
        datfile = twodF_DIR + '/data' + cred + creg + '_ntar.dat'
      else:
        if (iset < 9):
          cset = '00' + str(iset+1)
        else:
          cset = '0' + str(iset+1)
        datfile = twodF_DIR + '/rand' + cset + cred + creg + '_ntar.dat'
      print datfile

      f = open(datfile,'r')
      lines = f.readlines()[3:]

      for line in lines:
        fields = line.split()
        ras2df.append(float(fields[0]))
        dec2df.append(float(fields[1]))
        red2df.append(float(fields[2]))
        weifkp2df.append(float(fields[6]))
  f.close()
  ras2df,dec2df,red2df,weifkp2df = np.array(ras2df),np.array(dec2df),np.array(red2df),np.array(weifkp2df)
  n2df = len(ras2df)
  print n2df,'2dFLenS lenses'
  return ras2df,dec2df,red2df,weifkp2df,n2df

#============================================
# Find magnitudes/colours of closest source to each lens
def matchlenstosource(raslens,declens,rassource,decsource,grcolsource,ricolsource,rmagsource,rmin1,rmax1,dmin1,dmax1,rmin2,rmax2,dmin2,dmax2):
  print '\nFinding closest source to each lens...'
  separcmax = 1. # Matching separation in arcsec
  nlens = len(raslens)
  grcollens,ricollens,rmaglens,iphotlens = np.zeros(nlens),np.zeros(nlens),np.zeros(nlens),np.zeros(nlens,dtype='int')
  indexlens = np.arange(nlens)
  cut1 = (raslens > rmin1) & (raslens < rmax1) & (declens > dmin1) & (declens < dmax1)
  cut2 = ((raslens > rmin2+360.) | (raslens < rmax2)) & (declens > dmin2) & (declens < dmax2)
  cut = (cut1 | cut2)
  raslens1,declens1,indexlens1 = raslens[cut],declens[cut],indexlens[cut]
  print len(raslens1),'lenses in angular area'
  coosource = SkyCoord(rassource*u.deg,decsource*u.deg)
  coolens = SkyCoord(raslens1*u.deg,declens1*u.deg)
  indexsource,sep,d3d = coolens.match_to_catalog_sky(coosource)
  grcollens1,ricollens1,rmaglens1 = grcolsource[indexsource],ricolsource[indexsource],rmagsource[indexsource]
  cut = (sep.arcsec < separcmax)
  raslens1,declens1,grcollens1,ricollens1,rmaglens1,indexlens1 = raslens1[cut],declens1[cut], \
                                                                 grcollens1[cut],ricollens1[cut],\
                                                                 rmaglens1[cut],indexlens1[cut]
  nlens = len(raslens1)
  print nlens,'lenses matched within',separcmax,'arcsec'
  grcollens[indexlens1] = grcollens1
  ricollens[indexlens1] = ricollens1
  rmaglens[indexlens1] = rmaglens1
  iphotlens[indexlens1] = 1
  return grcollens,ricollens,rmaglens,iphotlens

#============================================
# Determine weights of catalogue to match magnitudes of reference using
# the KV450 DIR method
def calcmagweights(magscat,magsref,weiref):
  print '\nCalculating magnitude weights...'
  no_NN = 10
  ncat = magscat.shape[0]
  nref = magsref.shape[0]
# Build tree
  print '\nBuilding trees...'
  treecat = scipy.spatial.cKDTree(magscat,leafsize=100)
  treeref = scipy.spatial.cKDTree(magsref,leafsize=100)
# Nearest catalogue neighbours to each catalogue object
  neighbours_cat_of_cat = ( treecat.query(magscat,k=no_NN) )
  average_ref_weight = np.average(weiref)
  no_neighbours_ref_of_cat = np.zeros(ncat)
  neighbours_ref_of_cat = []
  weight_ref_of_cat = np.zeros(ncat)
  weicat = np.ones(ncat)
# Loop over each catalogue object
  for i in range(ncat):
# Indices of nearest reference neighbours to each catalogue object
    x = magscat[i,:]
    r = neighbours_cat_of_cat[0][i,no_NN-1]
    iref = treeref.query_ball_point(x,r)
    neighbours_ref_of_cat.append(iref)
    no_neighbours_ref_of_cat[i] = float(len(neighbours_ref_of_cat[i]))
    if (no_neighbours_ref_of_cat[i] > 0.):
      weight_ref_of_cat[i] = (np.average(weiref[neighbours_ref_of_cat[i]]))
      weicat[i] = (
                    (float(ncat)/float(nref)) *
                    (weight_ref_of_cat[i]/average_ref_weight) *
                    (no_neighbours_ref_of_cat[i]/float(no_NN))
                  )
  print len(no_neighbours_ref_of_cat[no_neighbours_ref_of_cat == 0.]),'catalogue objects with no neighbours'
  print 'Mean reference weight =',np.average(weiref)
  print 'Mean catalogue weight =',np.average(weicat)
  return weicat

#============================================
# Write out lens fits catalogue
def writelenscat(outfile,raslens,declens,redlens,weicomplens,weifkplens,iphotlens,weimaglens,grcollens,ricollens,rmaglens):
  print '\nWriting out lens catalogue...'
  print outfile
  col1 = fits.Column(name='ALPHA_J2000',format='D',array=raslens)
  col2 = fits.Column(name='DELTA_J2000',format='D',array=declens)
  col3 = fits.Column(name='Z',format='E',array=redlens)
  col4 = fits.Column(name='WEICOMP',format='E',array=weicomplens)
  col5 = fits.Column(name='WEIFKP',format='E',array=weifkplens)
  col6 = fits.Column(name='FLAGPHOT',format='J',array=iphotlens)
  col7 = fits.Column(name='WEIMAG',format='E',array=weimaglens)
  col8 = fits.Column(name='GRCOL',format='E',array=grcollens)
  col9 = fits.Column(name='RICOL',format='E',array=ricollens)
  col10 = fits.Column(name='RMAG',format='E',array=rmaglens)
  hdulist = fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10])
  hdulist.writeto(outfile)
  return

#============================================
# Read in lens fits catalogue
def readlenscat(infile):
  print '\nReading in lens catalogue...'
  print infile
  hdulist = fits.open(infile)
  table = hdulist[1].data
  raslens = table.field('ALPHA_J2000')
  declens = table.field('DELTA_J2000')
  redlens = table.field('Z')
  weicomplens = table.field('WEICOMP')
  weifkplens = table.field('WEIFKP')
  iphotlens = table.field('FLAGPHOT')
  weimaglens = table.field('WEIMAG')
  grcollens = table.field('GRCOL')
  ricollens = table.field('RICOL')
  rmaglens = table.field('RMAG')
  hdulist.close()
  nlens = len(raslens)
  print 'Read in',nlens,'lenses'
  return raslens,declens,redlens,weicomplens,weifkplens,iphotlens,weimaglens,grcollens,ricollens,rmaglens,nlens

#============================================
# Run test plots of the lens catalogues
def testcats(ired):
  opt = 4 # 1) (R.A., Dec.) overplot
          # 2) redshift overplot
          # 3) magnitude files
          # 4) weighted distributions
          # 5) weights
# Read in lens fits catalogues
#  stem = '/Users/cblake/Data/kids1000/lenscats/'
  stem = '/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/'
  if (ired == 1):
    cred = '_bz1'
    zmin,zmax = 0.2,0.5
  elif (ired == 2):
    cred = '_bz2'
    zmin,zmax = 0.4,0.6
  else:
    cred = '_bz3'
    zmin,zmax = 0.5,0.75
  infile = stem + 'boss_data_lenses' + cred + '.fits'
  rasbossdat,decbossdat,redbossdat,weicompbossdat,weifkpbossdat,iphotbossdat,weimagbossdat,grcolbossdat,ricolbossdat,rmagbossdat,nbossdat = readlenscat(infile)
  infile = stem + 'boss_random_lenses' + cred + '.fits'
  rasbossran,decbossran,redbossran,weicompbossran,weifkpbossran,iphotbossran,weimagbossran,grcolbossran,ricolbossran,rmagbossran,nbossran = readlenscat(infile)
  infile = stem + '2dflens_data_lenses' + cred + '.fits'
  ras2dfdat,dec2dfdat,red2dfdat,weicomp2dfdat,weifkp2dfdat,iphot2dfdat,weimag2dfdat,grcol2dfdat,ricol2dfdat,rmag2dfdat,n2dfdat = readlenscat(infile)
  infile = stem + '2dflens_random_lenses' + cred + '.fits'
  ras2dfran,dec2dfran,red2dfran,weicomp2dfran,weifkp2dfran,iphot2dfran,weimag2dfran,grcol2dfran,ricol2dfran,rmag2dfran,n2dfran = readlenscat(infile)
# Overplot data and random lenses
  if (opt == 1):
    ras1,dec1,lab1 = rasbossdat,decbossdat,'Data'
    ras2,dec2,lab2 = rasbossran,decbossran,'Randoms'
#    ras1,dec1,lab1 = ras2dfdat,dec2dfdat,'Data'
#    ras2,dec2,lab2 = ras2dfran,dec2dfran,'Randoms'
#    ras1,dec1,lab1 = rasbossdat,decbossdat,'BOSS'
#    ras2,dec2,lab2 = ras2dfdat,dec2dfdat,'2dFLenS'
#    ras1[ras1 > 180.] = ras1[ras1 > 180.] - 360.
#    ras2[ras2 > 180.] = ras2[ras2 > 180.] - 360.
#    rmin,rmax,dmin,dmax = 90.,270.,-15.,90.
#    rmin,rmax,dmin,dmax = -90.,90.,-90.,-25.
#    cut = (ras1 > rmin) & (ras1 < rmax) & (dec1 > dmin) & (dec1 < dmax)
#    ras1,dec1 = ras1[cut],dec1[cut]
#    cut = (ras2 > rmin) & (ras2 < rmax) & (dec2 > dmin) & (dec2 < dmax)
#    ras2,dec2 = ras2[cut],dec2[cut]
    fig = plt.figure()
    n1,n2 = len(ras1),len(ras2)
    nplot = min(10000,n1,n2)
    if (n1 > nplot):
      cut = np.random.choice(n1,nplot,replace=False)
    else:
      cut = np.full(n1,True,dtype=bool)
    plt.scatter(ras1[cut],dec1[cut],s=0.5,marker='o',color='black',alpha=0.25,label=lab1)
    if (n2 > nplot):
      cut = np.random.choice(n2,nplot,replace=False)
    else:
      cut = np.full(n2,True,dtype=bool)
    plt.scatter(ras2[cut],dec2[cut],s=0.5,marker='o',color='red',alpha=0.25,label=lab2)
    plt.xlabel('R.A. [deg]')
    plt.ylabel('Dec. [deg]')
    plt.legend()
    plt.show()
    sys.exit()
  elif (opt == 2):
    nz = 100
    red1,lab1 = redbossdat,'Data'
    red2,lab2 = redbossran,'Randoms'
#    red1,lab1 = red2dfdat,'Data'
#    red2,lab2 = red2dfran,'Randoms'
#    red1,lab1 = redbossdat,'BOSS'
#    red2,lab2 = red2dfdat,'2dFLenS'
    fig = plt.figure()
    hist1,zlims = np.histogram(red1,bins=nz,range=[zmin,zmax],normed=True)
    zcen = zlims[:-1] + 0.5*(zmax-zmin)/nz
    plt.plot(zcen,hist1,color='black',label=lab1)
    hist2,zlims = np.histogram(red2,bins=nz,range=[zmin,zmax],normed=True)
    plt.plot(zcen,hist2,color='red',label=lab2)
    plt.xlabel('z')
    plt.ylabel('p(z)')
    plt.legend()
    plt.show()
    sys.exit()
# Write out matched catalogues
  elif (opt == 3):
    if (ired == 1):
      fileboss = 'phot_bossz1.dat'
      file2df = 'phot_2dflz1.dat'
      outfile2df = 'weights_2dflz1.dat'
    elif (ired == 2):
      fileboss = 'phot_bossz2.dat'
      file2df = 'phot_2dflz2.dat'
      outfile2df = 'weights_2dflz2.dat'
    elif (ired == 3):
      fileboss = 'phot_bossz3.dat'
      file2df = 'phot_2dflz3.dat'
      outfile2df = 'weights_2dflz3.dat'
    print fileboss
    f = open(fileboss,'w')
    for i in range(nbossdat):
      if (iphotbossdat[i] > 0):
        f.write('{} {} {} {} {} {} {}'.format(rasbossdat[i],decbossdat[i],redbossdat[i],weicompbossdat[i],grcolbossdat[i],ricolbossdat[i],rmagbossdat[i]) + '\n')
    f.close()
    print file2df
    f = open(file2df,'w')
    for i in range(n2dfdat):
      if (iphot2dfdat[i] > 0):
        f.write('{} {} {:7.5f} {} {} {} {}'.format(ras2dfdat[i],dec2dfdat[i],red2dfdat[i],weicomp2dfdat[i],grcol2dfdat[i],ricol2dfdat[i],rmag2dfdat[i]) + '\n')
    f.close()
    print outfile2df
    f = open(outfile2df,'w')
    f.write('# R.A., Dec., redshift, weight\n')
    for i in range(n2dfdat):
      if (iphot2dfdat[i] > 0):
        f.write('{} {} {:7.5f} {}'.format(ras2dfdat[i],dec2dfdat[i],red2dfdat[i],weimag2dfdat[i]) + '\n')
    f.close()
  elif (opt == 4):
    iphotcat,redcat,grcolcat,ricolcat,rmagcat,weimagcat = iphot2dfdat,red2dfdat,grcol2dfdat,ricol2dfdat,rmag2dfdat,weimag2dfdat
    iphotref,redref,grcolref,ricolref,rmagref,weicompref = iphotbossdat,redbossdat,grcolbossdat,ricolbossdat,rmagbossdat,weicompbossdat
#    iphotcat,redcat,grcolcat,ricolcat,rmagcat,weimagcat = iphotbossdat,redbossdat,grcolbossdat,ricolbossdat,rmagbossdat,weimagbossdat
#    iphotref,redref,grcolref,ricolref,rmagref,weicompref = iphot2dfdat,red2dfdat,grcol2dfdat,ricol2dfdat,rmag2dfdat,weicomp2dfdat
    cutref = (iphotref > 0)
    cutcat = (iphotcat > 0)
    norm = np.sum(weicompref[cutref])/np.sum(weimagcat[cutcat])
    normed = False
    label1,label2,label3 = 'ref','cat','cat weighted'
    fig = plt.figure()
    nrow,ncol = 2,2
    sub = fig.add_subplot(nrow,ncol,1)
    xmin,xmax = 0.5,2.5
    sub.hist(grcolref[cutref],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='black',label=label1)
    sub.hist(grcolcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='red',label=label2)
    sub.hist(grcolcat[cutcat],weights=norm*weimagcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='blue',label=label3)
    sub.set_xlabel('g-r')
    sub.set_xlim(xmin,xmax)
    ymin,ymax = sub.get_ylim()
    sub.set_ylim(0.,1.1*ymax)
    plt.legend(prop={'size':10},loc=2)
    sub = fig.add_subplot(nrow,ncol,2)
    xmin,xmax = 0.,1.5
    sub.hist(ricolref[cutref],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='black',label=label1)
    sub.hist(ricolcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='red',label=label2)
    sub.hist(ricolcat[cutcat],weights=norm*weimagcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='blue',label=label3)
    sub.set_xlabel('r-i')
    sub.set_xlim(xmin,xmax)
    ymin,ymax = sub.get_ylim()
    sub.set_ylim(0.,1.1*ymax)
    sub = fig.add_subplot(nrow,ncol,3)
    xmin,xmax = 16.,23.
    sub.set_xlabel('r')
    sub.hist(rmagref[cutref],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='black',label=label1)
    sub.hist(rmagcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='red',label=label2)
    sub.hist(rmagcat[cutcat],weights=norm*weimagcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='blue',label=label3)
    sub.set_xlim(xmin,xmax)
    ymin,ymax = sub.get_ylim()
    sub.set_ylim(0.,1.1*ymax)
    sub = fig.add_subplot(nrow,ncol,4)
    xmin,xmax = zmin,zmax
    sub.hist(redref[cutref],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='black',label=label1)
    sub.hist(redcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='red',label=label2)
    sub.hist(redcat[cutcat],weights=norm*weimagcat[cutcat],bins=100,range=[xmin,xmax],histtype='step',normed=normed,facecolor='None',edgecolor='blue',label=label3)
    sub.set_xlabel('Redshift')
    sub.set_xlim(xmin,xmax)
    ymin,ymax = sub.get_ylim()
    sub.set_ylim(0.,1.1*ymax)
    fig.tight_layout()
    plt.show()
    sys.exit()
  elif (opt == 5):
    wmin,wmax,nw = -0.1,5.1,100
    wei1,lab1 = weicompbossdat,'BOSS completeness weight'
    wei2,lab2 = weifkpbossdat,'BOSS FKP weight'
    wei3,lab3 = weimagbossdat,'BOSS magnitude weight'
#    wei1,lab1 = weicomp2dfdat,'2dFLenS completeness weight'
#    wei2,lab2 = weifkp2dfdat,'2dFLenS FKP weight'
#    wei3,lab3 = weimag2dfdat,'2dFLenS magnitude weight'
#    wei1,lab1 = weicompbossran,'BOSS completeness weight'
#    wei2,lab2 = weifkpbossran,'BOSS FKP weight'
#    wei3,lab3 = weimagbossran,'BOSS magnitude weight'
#    wei1,lab1 = weicomp2dfran,'2dFLenS completeness weight'
#    wei2,lab2 = weifkp2dfran,'2dFLenS FKP weight'
#    wei3,lab3 = weimag2dfran,'2dFLenS magnitude weight'
    fig = plt.figure()
    hist1,lims = np.histogram(wei1,bins=nw,range=[wmin,wmax],normed=True)
    wcen = lims[:-1] + 0.5*(wmax-wmin)/nw
    plt.plot(wcen,hist1,color='black',label=lab1)
    hist2,lims = np.histogram(wei2,bins=nw,range=[wmin,wmax],normed=True)
    plt.plot(wcen,hist2,color='red',label=lab2)
    hist3,lims = np.histogram(wei3,bins=nw,range=[wmin,wmax],normed=True)
    plt.plot(wcen,hist3,color='blue',label=lab3)
    plt.xlabel('weight')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()
    sys.exit()
  return

#===================================
# We're now ready to run the script with command line options
# To do this in main we would need to define them as global

# Read in user input to set the location of the input/output and the desired bin
if len(sys.argv) <6: 
    #print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? lenscat.fits \
    #        randomcat.fits sourcecat.fits outfilename" % sys.argv[0]) 
  print "Usage: %s lens_bin KiDS_Location KiDS_Version BOSS_Location 2dFLenS_Location Out_Directory" % sys.argv[0] 
  print "Example python makelenscats.py 1 /disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/ rband_mask_BRIGHT_v3.cat \
          /disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/BOSS_original /disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/2dFLenS_original \
          /disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/"
  sys.exit(1)
else:
  ired = int(sys.argv[1]) 
  KiDS_DIR = sys.argv[2]
  KiDS_VER = sys.argv[3]
  BOSS_DIR = sys.argv[4]
  twodF_DIR = sys.argv[5]
  OUTDIR = sys.argv[6]

# Redshift bin for catalogues
# These are fixed by the Sanchez et al paper and so
# we hardwire these properties here
#  ired #(1) 0.2-0.5 (2) 0.5-0.75 (3) overlap - 0.4-0.6 - not used

# Generate the lens catalogues
makecats(ired)

# Run test plots of the lens catalogues
#  testcats(ired)
