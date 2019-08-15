# 05/03/2019, B. M. Giblin, PhD Student, Edinburgh
# Read in PSF residual quantities & calculate rho statistics

import numpy as np
from astropy.io import fits
import treecorr
import glob
import os
import sys

Calc_rhos = True         # Runs TreeCorr to calc rho1-5 if True.
PreReadCat = False       # If True, read in a pre-saved shear catalogue
Splitup_Fields = True    # If True, splits sky survey into Res*Res patches
Res = 7                  # Splits the Field into Res*Res pieces and calculates rho's for each.
                         # These are used to calculate Jackknife errors.

LFver = sys.argv[1]        # e.g. "321" or "309b" # Lensfit versions
NorS = sys.argv[2]         # "N" or "S" for North/South Fields
print("LFver is ",LFver)
print("NorS is ",NorS)
if LFver == "309b":
    expname = 'PSFRES_XI_svn_%s' %LFver
else:
    expname = 'PSFRES_XI_glab_%s' %LFver
    
#Fields = ["KIDS_129p0_0p5"]     # Set to this to read just one Field

# DECIDE WHICH FIELDS TO USE
#fname = '/home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/G9.txt'   # JUST USE G9
fname = '/home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/K1000_%s.txt' %NorS # USE ALL OF N OR S.
with open(fname) as f:
    Fields_tmp = f.readlines()
Fields=[]
for x in Fields_tmp:
    Fields.append(x.split()[1])

data_indir = "/disk09/KIDS/KIDSCOLLAB_V1.0.0/" # Contains all the KiDS Field directories
outdir = "/home/bengib/KiDS1000_NullTests/PSF_rho_stats/LFver%s" %LFver
if not os.path.exists(outdir):
    os.makedirs(outdir)
    os.makedirs(outdir+"/rho1")
    os.makedirs(outdir+"/rho2")
    os.makedirs(outdir+"/rho3")
    os.makedirs(outdir+"/rho4")
    os.makedirs(outdir+"/rho5")
    os.makedirs(outdir+"/Catalogues")


if PreReadCat:
    RA, Dec, e0PSF, e1PSF, delta_e0PSF, delta_e1PSF, TPSF, delta_TPSF, Xpos, Ypos = np.load('%s/Catalogues/PSF_Data_%s.npy'%(outdir,NorS)).transpose()
else:
    Field_ID = np.empty([]) # ID of the field in file 'fname'
                            # where the PSF data came from
    # RA, dec
    RA = np.empty([])
    Dec = np.empty([])
    # PSF ellip. & residuals
    e0PSF = np.empty([])
    e1PSF = np.empty([])
    delta_e0PSF = np.empty([])
    delta_e1PSF = np.empty([])
    # PSF size
    TPSF = np.empty([])
    delta_TPSF = np.empty([])
    Xpos = np.empty([])       # These guys aren't used for rho stat cal
    Ypos = np.empty([])       # But they need to be saved to plot PSF residual across chip over on eday.
    i=0
    for f in Fields:
        print("Reading in field %s of %s"%(i,len(Fields)))
        exposures = glob.glob('%s/%s/checkplots/%s/*_PSFres.cat'%(data_indir,f,expname))
        for e in exposures:
            Field_ID = np.append(Field_ID,i)
            fitsfile = fits.open(e)
            tmp_RA = fitsfile[1].data['ALPHA_J2000']
            tmp_Dec = fitsfile[1].data['DELTA_J2000']
            tmp_Xpos = fitsfile[1].data['Xpos']
            tmp_Ypos = fitsfile[1].data['Ypos']
            RA = np.append(RA, tmp_RA)
            Dec = np.append(Dec, tmp_Dec)
            Xpos = np.append(Xpos, tmp_Xpos)
            Ypos = np.append(Ypos, tmp_Ypos)
            # Moments
            tmp_I_xx = fitsfile[1].data['moments_data0']
            tmp_I_yy = fitsfile[1].data['moments_data1']    
            tmp_I_xy = fitsfile[1].data['moments_data2']    
            tmp_I_xx_mod = fitsfile[1].data['moments_model0']
            tmp_I_yy_mod = fitsfile[1].data['moments_model1'] 
            tmp_I_xy_mod = fitsfile[1].data['moments_model2'] 
            tmp_TPSF = tmp_I_xx + tmp_I_yy
            tmp_TPSF_mod = tmp_I_xx_mod + tmp_I_yy_mod
            # PSF sizes
            TPSF = np.append(TPSF, tmp_TPSF)
            delta_TPSF = np.append(delta_TPSF, tmp_TPSF - tmp_TPSF_mod)
            
            # PSF ellipticities
            # First for data
            tmp_e0PSF = fitsfile[1].data['psfe_data0']
            tmp_e1PSF = fitsfile[1].data['psfe_data1']
            e0PSF = np.append(e0PSF, tmp_e0PSF)
            e1PSF = np.append(e1PSF, tmp_e1PSF)
            # Now the residuals (data - model)
            tmp_e0PSF_mod = fitsfile[1].data['psfe_model0']
            tmp_e1PSF_mod = fitsfile[1].data['psfe_model1']
            tmp_delta_e0PSF = tmp_e0PSF - tmp_e0PSF_mod
            tmp_delta_e1PSF = tmp_e1PSF - tmp_e1PSF_mod
            delta_e0PSF = np.append(delta_e0PSF, tmp_delta_e0PSF)
            delta_e1PSF = np.append(delta_e1PSF, tmp_delta_e1PSF)
        i+=1
    # There's something weird with element 0: RA=745... Delete it.
    RA = np.delete(RA, 0)
    Dec = np.delete(Dec, 0)
    Xpos = np.delete(Xpos, 0)
    Ypos = np.delete(Ypos, 0)
    e0PSF = np.delete(e0PSF, 0)
    e1PSF = np.delete(e1PSF, 0)
    delta_e0PSF = np.delete(delta_e0PSF, 0)
    delta_e1PSF = np.delete(delta_e1PSF, 0)
    TPSF = np.delete(TPSF, 0)
    delta_TPSF = np.delete(delta_TPSF, 0)

    colstack = np.column_stack((RA, Dec,
                                e0PSF, e1PSF, delta_e0PSF, delta_e1PSF, TPSF, delta_TPSF, Xpos,Ypos))
    np.save('%s/Catalogues/PSF_Data_%s.npy'%(outdir,NorS), colstack)


ThBins=9
min_sep=0.5   # arcmin
max_sep=300.  # arcmin                                                                                        
bin_slop=0.1/(np.log(max_sep/min_sep)/float(ThBins))
metric='Arc'

def Run_TreeCorr(ra,dec,y1,y2,z1,z2):
    cat1 = treecorr.Catalog(ra=ra, dec=dec,
                ra_units='degrees', dec_units='degrees',
		g1=y1, g2=y2)
    cat2 = treecorr.Catalog(ra=ra, dec=dec,
                ra_units='degrees', dec_units='degrees',
                g1=z1, g2=z2)    
    gg = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep,
                bin_slop=bin_slop, nbins=ThBins,
                metric=metric, sep_units='arcmin')
    gg.process(cat1,cat2)
    return gg.meanr, gg.xip, gg.xim




# This one we don't use - returns elements INSIDE (ra,dec) patch
def Select_Patch(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the elements in Q corresponding to INSIDE the (ra,dec) range 
    idx_ra = np.where(np.logical_and(ra<rhi, ra>rlo))[0]
    idx_dec = np.where(np.logical_and(dec<dhi, dec>dlo))[0]
    idx = np.intersect1d(idx_ra, idx_dec) 
    return Q[idx]


def Split_Fields(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the elements in Q corresponding to OUTSIDE the (ra,dec) range (Jackknife)
    idx_ra = np.append( np.where(ra>rhi)[0], np.where(ra<rlo)[0] )
    idx_dec = np.append( np.where(dec>dhi)[0], np.where(dec<dlo)[0] )
    idx = np.unique( np.append(idx_ra, idx_dec) ) # delete the things that appear more than once
    return Q[idx]



if Calc_rhos:

    if Splitup_Fields:
        print("Splitting up %s Field into %s*%s patches" %(NorS,Res,Res))
        if True in (RA>300) and True in (RA<10):
            # RA's cross RA=0, must be more canny in defining RA_binning
            RA[ RA>300 ] -= 360 # Not a general solution; works for KiDS-S as RA[RA>300].max() = 328  
        
        RA_coarse = np.linspace(RA.min(), RA.max(), Res+1)
        Dec_coarse = np.linspace(Dec.min(), Dec.max(), Res+1)
        for i in range(Res*Res):
            print("On patch %s of %s" %(i,Res*Res))
            ridx = i % Res
            didx = int( i / Res )
            rlo = RA_coarse[ridx]
            rhi = RA_coarse[ridx+1]
            dlo = Dec_coarse[didx]
            dhi = Dec_coarse[didx+1]
            tmp_RA = Split_Fields(RA, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_Dec = Split_Fields(Dec, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_e0PSF = Split_Fields(e0PSF, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_e1PSF = Split_Fields(e1PSF, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_delta_e0PSF = Split_Fields(delta_e0PSF, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_delta_e1PSF = Split_Fields(delta_e1PSF, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_TPSF = Split_Fields(TPSF, RA, Dec, rlo, rhi, dlo, dhi)
            tmp_delta_TPSF = Split_Fields(delta_TPSF, RA, Dec, rlo, rhi, dlo, dhi)
            if len(tmp_RA) > 100: # some patches will be empty
                TScaler = (tmp_delta_TPSF/tmp_TPSF)
                print("Calculating rho statistics")
                print("On rho1")
                meanr1, rhop1, rhom1 = Run_TreeCorr( tmp_RA, tmp_Dec, tmp_delta_e0PSF, tmp_delta_e1PSF, tmp_delta_e0PSF, tmp_delta_e1PSF )
                print("On rho2")
                meanr2, rhop2, rhom2 = Run_TreeCorr( tmp_RA, tmp_Dec, tmp_e0PSF, tmp_e1PSF, tmp_delta_e0PSF, tmp_delta_e1PSF )
                print("On rho3")
                meanr3, rhop3, rhom3 = Run_TreeCorr( tmp_RA, tmp_Dec, tmp_e0PSF*TScaler, tmp_e1PSF*TScaler, tmp_e0PSF*TScaler, tmp_e1PSF*TScaler )
                print("On rho4")
                meanr4, rhop4, rhom4 = Run_TreeCorr( tmp_RA, tmp_Dec, tmp_delta_e0PSF, tmp_delta_e1PSF, tmp_e0PSF*TScaler, tmp_e1PSF*TScaler )
                print("On rho5")
                meanr5, rhop5, rhom5 = Run_TreeCorr( tmp_RA, tmp_Dec, tmp_e0PSF, tmp_e1PSF, tmp_e0PSF*TScaler, tmp_e1PSF*TScaler )
                np.savetxt('%s/rho1/rho1_KiDS_%s_%sof%sx%s.dat'%(outdir,NorS,i,Res,Res), np.c_[meanr1, rhop1, rhom1],
                           header='# mean-theta [arcmin], rho1_p, rho1_m')
                np.savetxt('%s/rho2/rho2_KiDS_%s_%sof%sx%s.dat'%(outdir,NorS,i,Res,Res), np.c_[meanr2, rhop2, rhom2],
                           header='# mean-theta [arcmin], rho2_p, rho2_m')
                np.savetxt('%s/rho3/rho3_KiDS_%s_%sof%sx%s.dat'%(outdir,NorS,i,Res,Res), np.c_[meanr3, rhop3, rhom3],
                           header='# mean-theta [arcmin], rho3_p, rho3_m')
                np.savetxt('%s/rho4/rho4_KiDS_%s_%sof%sx%s.dat'%(outdir,NorS,i,Res,Res), np.c_[meanr4, rhop4, rhom4],
                           header='# mean-theta [arcmin], rho4_p, rho4_m')
                np.savetxt('%s/rho5/rho5_KiDS_%s_%sof%sx%s.dat'%(outdir,NorS,i,Res,Res), np.c_[meanr5, rhop5, rhom5],
                           header='# mean-theta [arcmin], rho5_p, rho5_m')

                
                        
    else:
        print("Calculating rho statistics")
        print("On rho1")
        meanr1, rhop1, rhom1 = Run_TreeCorr( RA, Dec, delta_e0PSF, delta_e1PSF, delta_e0PSF, delta_e1PSF )
        print("On rho2")
        meanr2, rhop2, rhom2 = Run_TreeCorr( RA, Dec, e0PSF, e1PSF, delta_e0PSF, delta_e1PSF )
        print("On rho3")
        meanr3, rhop3, rhom3 = Run_TreeCorr( RA, Dec, e0PSF*(delta_TPSF/TPSF), e1PSF*(delta_TPSF/TPSF), e0PSF*(delta_TPSF/TPSF), e1PSF*(delta_TPSF/TPSF) )
        print("On rho4")
        meanr4, rhop4, rhom4 = Run_TreeCorr( RA, Dec, delta_e0PSF, delta_e1PSF, e0PSF*(delta_TPSF/TPSF), e1PSF*(delta_TPSF/TPSF) )
        print("On rho5")
        meanr5, rhop5, rhom5 = Run_TreeCorr( RA, Dec, e0PSF, e1PSF, e0PSF*(delta_TPSF/TPSF), e1PSF*(delta_TPSF/TPSF) )

        np.savetxt('%s/rho1/rho1_KiDS_%s.dat'%(outdir,NorS), np.c_[meanr1, rhop1, rhom1],
                   header='# mean-theta [arcmin], rho1_p, rho1_m')
        np.savetxt('%s/rho2/rho2_KiDS_%s.dat'%(outdir,NorS), np.c_[meanr2, rhop2, rhom2],
           header='# mean-theta [arcmin], rho2_p, rho2_m')
        np.savetxt('%s/rho3/rho3_KiDS_%s.dat'%(outdir,NorS), np.c_[meanr3, rhop3, rhom3],
           header='# mean-theta [arcmin], rho3_p, rho3_m')
        np.savetxt('%s/rho4/rho4_KiDS_%s.dat'%(outdir,NorS), np.c_[meanr4, rhop4, rhom4],
           header='# mean-theta [arcmin], rho4_p, rho4_m')
        np.savetxt('%s/rho5/rho5_KiDS_%s.dat'%(outdir,NorS), np.c_[meanr5, rhop5, rhom5],
           header='# mean-theta [arcmin], rho5_p, rho5_m')
