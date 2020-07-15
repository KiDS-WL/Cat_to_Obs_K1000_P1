#Script from Hendrik - hacked by Catherine
# Hacked again by B. Giblin in the winter of 2019/20 and built
# better, faster, stronger. We have the technology.

# The code now works on MICE mocks, K1000/BOSS/GAMA data,
# and the gamma_t model includes a new parameter to
# fit for the effects of magnification.
# Also you can toggle Single_Bin True/False to perform
# SRT for a single lens+source bin (specified on the command line).

# Now we use treecorr and Mandelbaum estimator
# to do the random and boost correction
# I also removed the scale selection part
# and we just use BOSS, so it's easier than
# the previous 2 spec-catalogue

import matplotlib.pyplot as plt
import numpy as np
import sys
from math import sqrt
from scipy.optimize import curve_fit, minimize
from scipy import stats
import string
import matplotlib.ticker as ticker
import string
import numpy.ma as ma
from astropy.io import ascii
from scipy.optimize import curve_fit
from astropy.io import fits

plt.rc('font', size=10)
# ----- Load Input ----- #
DFLAG = ''                     # Data Flag. Set to '' for Fiducial MICE. '_Octant' for MICE Octant.
Include_mCov = True            # Include the uncertainty due to the m-correction
Include_mBias = False           # If True, read in and bias the gamma_t by *(1 + m + f_mBias * m_err)
f_mBias = 5                    # where f_mBias (no. of sigmas to bias by) is defined below.

Include_Hartlap = False        # Hartlap correction
Include_Magnification = False   # If True, include extra param in gamma_t model: strength of magnifcation effect on gt.

Include_IA = True              # If True, read in the kcap prediction for the IA-only <g_t> and inflate the cov diag by
f_IA = 0                       # +(f_IA* gt_IA)^2 ; f_IA is our uncert. on IA amplitude.
A_IA = 1                       # A_IA*gt_IA is added to the fitted model if Include_IA is True.

Single_Bin = False             # If True, fit for only a single lens-source bin, specified by user on the command line.
                               # Else fit for all bins simultaneously.
                               
nofz_shift="_ModnofzUpDown"       # Only for K1000: use the Dls/Ds values for the nofz which has been
                               # shifted up ('_nofzUp'), down ('_nofzDown') by +/-(delta_z+delta_z_err)
                               # OR 5sig shift-up ('_nofzUp5sig'), down ('_nofzDown5sig') by (+/- 5*delta_z_err)
                               # For no shift, set to ''.
                               # Finally, to include the uncert. on the nofz's in the SRT modelling,
                               # set nofz_shift to "_ModnofzUpDown"
                               
from Get_Input import Get_Input
paramfile = sys.argv[1]   # e.g. params_KiDSdata.dat
GI = Get_Input(paramfile)
SOURCE_TYPE = GI.Source_Type()
LENS_TYPE = GI.Lens_Type()
RANDOM_TYPE = GI.Random_Type()

if Single_Bin == True:
    print("Single_Bin is True, so fitting one amplitude parameter to the specified lens and source bin.")
    print("These must be specified, lens-bin(0-4), source-bin(0-4), on the command line.")
    spec_bins = [ int(sys.argv[2]) ] # spec bin number (0-4)
    tomo_bins = [ int(sys.argv[3]) ] # tomo bin number (0-4)
    ntomo=1                     # Number of photo-z bins
    nspecz=1                    # "..." of spec-z bins

    # Note the tomo & source bin in the output save-file
    save_tag = '_t%ss%s' %(tomo_bins[0],spec_bins[0])

else:
    print("Performing SRT for all lens and source bins simultaneously. ",
          "Will fit an amplitude parameter per lens and source bin.")
    spec_bins=range(5)
    tomo_bins=range(5)
    ntomo=len(tomo_bins)              # Number of photo-z bins
    nspecz=len(spec_bins)             # "..." of spec-z bins
    if ntomo==5 and nspecz==5:
        save_tag = ''                 # Blank save_tag means SRT was conducted on all lens & source bins
    else:
        save_tag = '_%stbins-%ssbins' %(ntomo,nspecz) # Using more than a single t&s bin,
                                                     # but not all 5x5
    
numz_tag=5               # Part of 'GT_' filenames; GT_6Z_source_Y_${numz_tag}Z_lens_X.asc
ntheta = 4
theta_min = 2.
theta_max = 30.

measurements = {}

thetas_list = []  
gt_list = []      # gamma_t
gx_list = []      # gamma_x
gterr_list = []

# The following two are used to construct the analytical covariance
# and are only read in if Cov_Method="Analytical"
weight_sqrd_list = []  # the weight returned from the lens-source correlation with squared weights
ran_weight_list = []

tomo_list = []
speczbin_list = []
Dls_over_Ds_list = []
Dls_over_Ds_list_UP = []   # These two lists only get used if you are reading in the Dls/Ds for BOTH the
Dls_over_Ds_list_DOWN = [] # nofz shift up/down, and including this extra uncert. in the fitting
                           # (i.e. only if nofz_shift is set to "_ModnofzUpDown")


Cov_Method = "Spin"   # The method for calculating the gamma_t realisations for use in covariance estimation
                       # "Spin" many spin realisations of the source ellipticities (ie - shape noise only)
                       # "Patch" using other MICE realisations (1/8 of the sky)
                       # divided them into patches and calcuted the gamma_t from each patch.
                       # "Analytical" means construct and use an analytical covariance matrix.
                       
nPatch = 16             # If Cov_Method is Patch, the MICE octant is split into nPatch RA
                       # and nPatch Dec slices (nPatch^2 patches in total). gamma_t is
                       # calculated from each patch.   


INDIR='Output/SOURCE-%s_LENS-%s' %(SOURCE_TYPE, LENS_TYPE)
DlsDIR='Dls_over_Ds_data/SOURCE-%s_LENS-%s' %(SOURCE_TYPE, LENS_TYPE)
if "MICE2" in SOURCE_TYPE:
    # Additional identifiers when working with MICE
    # True/Estimated P(z) and Magnification on/off
    Mag_OnOff = GI.Mag_OnOff()
    Pz = GI.Pz_TrueEstimated()
    SN = GI.SN()
    INDIR += '_Pz%s_SN%s_mag%s' %(Pz,SN,Mag_OnOff)
    DlsDIR += '_Pz%s_mag%s%s' %(Pz,Mag_OnOff,DFLAG)
elif "K1000" in SOURCE_TYPE:
    Blind = GI.Blind()
    SOMFLAGNAME = GI.SOMFLAGNAME()
    OL_Tag = GI.OL_Tag() # Whether a high-z outlier is included in n(z)'s
                         # that went into the Dls/Ds ratios
    INDIR += '_Blind%s_SOM%s' %(Blind,SOMFLAGNAME)
    DlsDIR += '_Blind%s_SOM%s%s' %(Blind,SOMFLAGNAME, OL_Tag)
    save_tag += OL_Tag

if Include_mBias:
    # Read in the m-bias and error per tomo_bin
    mBias, mBias_err = np.loadtxt('mbias_perbin_mean_err.asc', usecols=(0,1), unpack=True)
    #mBias_err[1] *= -1
    #mBias_err[3] *= -1
else:
    mBias = np.zeros(tomo_bins[-1]+1)
    mBias_err = np.zeros_like( mBias )
    
nspin=500
if Cov_Method == "Spin" or Cov_Method == "Analytical":
    ncycle = nspin
    OUTDIR = INDIR + "/SPIN/"
    Area_Scale = 1. #341./777. # Approx. scaling KV450-->K1000     #1.
elif Cov_Method == "Patch":
    ncycle = nPatch*nPatch
    OUTDIR = INDIR + "/PATCH/"
    Area_Scale = 4*np.pi *(180./np.pi)**2. / 343.
else:
    print("Cov_Method must be set to Spin or Patch. Currently it's set to %s. EXITING." %Cov_Method)
    sys.exit()

for tomobin in tomo_bins:
    for speczbin in spec_bins:
        print("Reading tomo bin %s, spec bin %s" %(tomobin+1, speczbin+1))
        gtfile='%s/GT_6Z_source_%s_%sZ_lens_%s%s.asc' %(INDIR,(tomobin+1), numz_tag,(speczbin+1), DFLAG)
        gtdat=ascii.read(gtfile)

        measurements['thetas_'+str(speczbin)+"_"+str(tomobin)] = gtdat['meanr']
        measurements['gt_'+str(speczbin)+"_"+str(tomobin)]     = gtdat['gamT']*( 1+(mBias+f_mBias*mBias_err)[tomobin] )
        measurements['gx_'+str(speczbin)+"_"+str(tomobin)]     = gtdat['gamX']*( 1+(mBias+f_mBias*mBias_err)[tomobin] )
        measurements['gterr_'+str(speczbin)+"_"+str(tomobin)]  = gtdat['sigma']*( 1+(mBias+f_mBias*mBias_err)[tomobin] )

        thetas_list.append(measurements['thetas_'+str(speczbin)+"_"+str(tomobin)])
        gt_list.append(measurements['gt_'+str(speczbin)+"_"+str(tomobin)])
        gx_list.append(measurements['gx_'+str(speczbin)+"_"+str(tomobin)])
        gterr_list.append(measurements['gterr_'+str(speczbin)+"_"+str(tomobin)])
        
        if Cov_Method == "Analytical":
            measurements['ranweight_'+str(speczbin)+"_"+str(tomobin)] = gtdat['ranweight']
            measurements['weight_sqrd_'+str(speczbin)+"_"+str(tomobin)] = gtdat['weight_sqrd']
            ran_weight_list.append(measurements['ranweight_'+str(speczbin)+"_"+str(tomobin)]) 
            weight_sqrd_list.append(measurements['weight_sqrd_'+str(speczbin)+"_"+str(tomobin)])
        
        tomo_list.append(tomobin*np.ones((ntheta),dtype=np.int16))
        speczbin_list.append(speczbin*np.ones((ntheta),dtype=np.int16))

        # Read in the Dls_over_Ds data created with Dls_over_Ds.py
        # Either read in ONLY the nofz shift up/down/no-shift...
        # OR read in all, and incorporate the uncertainty on mean-z into the fitting.
        Dls_over_Ds_file = '%s/Dls_over_Ds_DIR_6Z_source_%s_%sZ_lens_%s' %(DlsDIR, (tomobin+1), numz_tag,(speczbin+1))
        if nofz_shift == "_ModnofzUpDown":
            #Read in both nofz shift up and down (later use them to inflate the errors on the covariance):
            Dls_over_Ds_tmp_up = np.loadtxt(Dls_over_Ds_file+'_nofzUp.asc')
            Dls_over_Ds_tmp_down = np.loadtxt(Dls_over_Ds_file+'_nofzDown.asc')
            Dls_over_Ds_tmp = np.loadtxt(Dls_over_Ds_file+'.asc')
            
            measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)+"_UP"] = np.repeat(Dls_over_Ds_tmp_up, ntheta)
            measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)+"_DOWN"] = np.repeat(Dls_over_Ds_tmp_down, ntheta)
            measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)] = np.repeat(Dls_over_Ds_tmp, ntheta)

            Dls_over_Ds_list_UP.append(measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)+"_UP"])
            Dls_over_Ds_list_DOWN.append(measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)+"_DOWN"])
            Dls_over_Ds_list.append(measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)])
            
        else:
            Dls_over_Ds_tmp = np.loadtxt(Dls_over_Ds_file+'%s.asc' %nofz_shift)
            measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)] = np.repeat(Dls_over_Ds_tmp,ntheta)
            Dls_over_Ds_list.append(measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)])
        
thetas = np.hstack(thetas_list)
tomo = np.hstack(tomo_list)
speczbin = np.hstack(speczbin_list)
gt = np.hstack(gt_list) 
gx = np.hstack(gx_list) 
gterr = np.hstack(gterr_list)
Dls_over_Ds = np.hstack(Dls_over_Ds_list)
if nofz_shift == "_ModnofzUpDown":
    Dls_over_Ds_UP = np.hstack(Dls_over_Ds_list_UP)
    Dls_over_Ds_DOWN = np.hstack(Dls_over_Ds_list_DOWN)
    delta_A = gt * Dls_over_Ds * (1./Dls_over_Ds_UP - 1./Dls_over_Ds_DOWN) /2.0
    
nmatrix = nspecz*ntomo*ntheta
gtlens=np.zeros([ntheta,ncycle])
diag=np.zeros(nmatrix)
covdiag=np.zeros([nmatrix,nmatrix])
#INDIRcov='Output/SOURCE-MICE2_KV450_LENS-MICE2_BOSS_PzEstimated_SNTrue_magoff/'
INDIRcov=INDIR
for tomobin in tomo_bins:
    for speczbin in spec_bins:
        for icycle in range(ncycle):
            if Cov_Method == "Spin" or Cov_Method == "Analytical":
                gtfile='%s/SPIN/GT_SPIN_%s_6Z_source_%s_%sZ_lens_%s.asc' %(INDIRcov,icycle,
                                                                           tomobin+1,numz_tag,speczbin+1)
            elif Cov_Method == "Patch":
                gtfile='%s/PATCH/GT_PATCH_%sof%s_6Z_source_%s_%sZ_lens_%s.asc' %(INDIRcov,icycle,nPatch*nPatch,
                                                                                 tomobin+1,num_ztag,speczbin+1)
                
            gtcycledat=ascii.read(gtfile)
            gtlens[:,icycle]=gtcycledat['gamT'] 
        if speczbin==spec_bins[0] and tomobin==tomo_bins[0]:
            if Single_Bin:
                # Then we're only reading in one bin,
                gtall = np.copy(gtlens)
            else:
                # There's more bins to read in
                gtprev = np.copy(gtlens)
        else:    
            gtall = np.vstack((gtprev,gtlens))
            gtprev= np.copy(gtall)
            print(len(gtall))
cov_gt=np.cov(gtall)


if Cov_Method == "Analytical":
    # Construct the analytical covariance matrix
    weight_sqrd = np.hstack(weight_sqrd_list)
    ran_weight = np.hstack(ran_weight_list)

    # Need to read in the weights for the randoms and the lenses:
    N_oversample = np.zeros( numz_tag )
    for speczbin in spec_bins:
        lenscatname='LENSCATS/%s%s/lens_cat_%sZ_%s.fits' %(LENS_TYPE,DFLAG, numz_tag,speczbin+1)
        rancatname='LENSCATS/%s%s/lens_cat_%sZ_%s.fits' %(RANDOM_TYPE,DFLAG, numz_tag,speczbin+1)
        f_l = fits.open(lenscatname)
        w_l = f_l[1].data['WEICOMP']
        f_l.close()
        f_r = fits.open(rancatname)
        w_r = f_r[1].data['WEICOMP']
        f_r.close()
        N_oversample[speczbin] =np.sum(w_r)/np.sum(w_l)

    # Now calculate npairs_weighted for each tomo and spec bin
    sigma_e = np.loadtxt('sigma_e_AllBlinds.asc', usecols=(0,), unpack=True) # col (0,1,2) = blind (A,B,C)
    npairs_weighted = np.zeros_like( ran_weight )
    cov_gt_a = np.zeros([ len(ran_weight), len(ran_weight) ])
    for i in range(ntomo): #tomobin in tomo_bins:
        for j in range(nspecz): #speczbin in spec_bins:
            # Locate the elements of ran_weight and weight_sqrd to extract.
            minidx = i*nspecz*ntheta+ j*ntheta
            maxidx = minidx + ntheta
            npairs_weighted[ minidx:maxidx ] = ran_weight[minidx:maxidx]**2 / ( N_oversample[spec_bins[j]]**2 *weight_sqrd[minidx:maxidx] )
            # Fill the diagonal elements of the covariance matrix
            for t in range(minidx, maxidx):
                cov_gt_a[t,t] = sigma_e[spec_bins[j]]**2 / npairs_weighted[t]

    # Replace the spin or patch covariance with the analytical one.
    cov_gt = np.copy ( cov_gt_a )
            
    

if Include_mCov:
    # Here we include the uncertainty due to the m-correction
    sigma_m = [0.019, 0.020, 0.017, 0.012, 0.010] # m-uncert. per source bin.
    cov_m = np.zeros_like( cov_gt )
    
    # Arguably the gt in the eqn for m-cov should be theoretical, not noisy data.
    # So we can read in the fit we obtained when excluding m-corr. But this seems
    # to make fit slightly worse, so commented this out.
    #gt_theory = np.loadtxt('gt_theory.dat')
    for a in range( cov_m.shape[0] ):
        for b in range( cov_m.shape[1] ):
            # Find which tomo & spec bin element (a,b) corresponds to:

            # For first gamma_1
            i = int( a/(nspecz*ntheta) )   # tomo bin, changes every nspecz*ntheta elements
            #j = int( a % (ntomo*ntheta))  # spec bin, not necessary to calculate.
            
            # For second gamma_t....:
            k = int( b/(nspecz*ntheta) )
            #l = int( b % (ntomo*ntheta))  

            # The gt's here could be replaced by gt_theory if desired (doesnt seem to help at present).
            cov_m[a,b] = gt[a] * gt[b] * (sigma_m[i]*sigma_m[k])
else:
    cov_m = np.zeros_like(cov_gt)


cov = (cov_gt + cov_m) * Area_Scale
print(len(cov))


if Include_IA:
    IA_DIR = '/home/bengib/kcap_NewInst/kcap/examples/GGL_IA/output/output_%sx%s_%s%s/gt_binned_ia_only'%(SOURCE_TYPE,
                                                                                                          LENS_TYPE.split('_')[0],
                                                                                                          Blind,OL_Tag)
    gt_IA = np.zeros([ ntomo*nspecz, ntheta ])
    k=0
    for tomobin in tomo_bins:
        for speczbin in spec_bins:
            gt_IA[k,:] = np.loadtxt('%s/bin_%s_%s.txt' %(IA_DIR,speczbin+1,tomobin+1))
            k+=1
    gt_IA = gt_IA.flatten()
else:
    gt_IA = np.zeros_like(gt)        

# Inflate the covariance diagonal by the uncertainty introduced by the IA:
# This does nothing if Include_IA is False.
delta_IA = f_IA * gt_IA
    
for i in range(nmatrix):
    cov[i,i] += delta_IA[i]**2. * Area_Scale
    
    # If incl. the uncert. due to nofz shifts, inflate the diag of the covariance
    if nofz_shift == "_ModnofzUpDown":
        cov[i,i] += delta_A[i]**2 * Area_Scale
        
    diag[i]=np.sqrt(cov[i,i])
    covdiag[i,i]=cov[i,i]
        

    
Corr=np.corrcoef(gtall)
plt.imshow(Corr, interpolation='None')
plt.colorbar()
plt.axis('off')
#plt.show()
#plt.savefig(OUTDIR+'/%sx%s_Shear_ratio_correlation_matrix%s.png'%(SOURCE_TYPE,LENS_TYPE,save_tag))

if Single_Bin == False:
    # Save the full covariance matrix, so it can easily be read in by
    # Compare_BFParams_And_Models.py and errors easily extracted
    np.save(OUTDIR+'/GTCovMat_%sx%s_6Z_source_5Z_lens_mCov%s%s'%(SOURCE_TYPE,LENS_TYPE,Include_mCov,save_tag), cov)

#cov_inv=np.linalg.inv(cov)
# Catherine reckons these two lines are a more stable way of inverting
# the cov than line above in cases where it's close to singular.
# For 500 spin realisations & 5 lens bins, p-value is insensitive to which method we use.
U,s,V=np.linalg.svd(cov)
cov_inv = np.transpose(V) @ np.diag(1./s) @ np.transpose(U)
#cov_inv=np.linalg.inv(covdiag)

Hartlap = (ncycle - ntomo * nspecz * ntheta - 2.) / (ncycle - 1.)
if Include_Hartlap:
    cov_inv *= Hartlap

#############################

def func(params):
    amplitude = params[:(num_amplitudes)]
    amplitude = np.hstack((amplitude,)*ntomo) # repeat array ntomo times
    model = Dls_over_Ds * amplitude
    if Single_Bin:
        # If fitting only a single s&l bin, the amplitudes don't change with theta,
        # so we need to manually include the 1/theta dependence in the model.
        model *= (1./thetas)
    
    if Include_Magnification:
        # If including magnification there's an additional nspecz alpha free params.
        alpha = params[(num_amplitudes):]
        alpha = np.repeat(alpha, ntheta)  # repeat each element ntheta times (no. points per lens bin)
        alpha = np.hstack((alpha,)*ntomo) # then create ntomo copies of this array
        model = model + 2.*(alpha-1.) * Magnif_Shape

    if Include_IA:
        model += A_IA * gt_IA
        
    return model


######## including covariance ###

def chi2(params, data, cov_inv):
    if Include_Magnification:
        #params[num_amplitudes:] = np.ones(5) # Manually kill the magnification.
        
        # Apply a prior on the alpha parameters
        if params[num_amplitudes:].min()<0. or params[num_amplitudes:].max()>5.:
        #    print("min params is ", params[(ntheta*nspecz):].min(),
        #          " max params is ", params[(ntheta*nspecz):].max(),
        #          " returning chi2 of infinity.")
            return np.inf
    return np.dot(data-func(params),np.dot(cov_inv,data-func(params)))


# The gamma_t model depends on if we are running the SRT with all source & lens bins together
# Or just one source+lens bin.                                                                                                  
if Single_Bin:
    # There is only one free amplitude parameter.
    num_amplitudes = 1
else:
    # There are ntheta*nspecz amplitude free params                                                                             
    # (one per lens and theta bin).                                                                                             
    num_amplitudes = ntheta*nspecz
    
if Include_Magnification:
    # Read in the model for magnifcation shape, given the K1000/BOSS nofz's
    # (NB: this shape won't be a good fit for other nofz's)
    # There's also nspecz more free params in the model (one alpha per lens bin)
    Magnif_Shape = np.load('/home/bengib/kcap_NewInst/kcap/examples/output_magnification_alpha1.0/SRTparam_Bij.npy').flatten()
    
    if Single_Bin:
        # Extract just the elements of Magnif_Shape corresponding to the bins specified by the user
        idx_mag = tomo_bins[0]*ntheta*5 + spec_bins[0]*ntheta # 1st element in the specified bin
        Magnif_Shape = Magnif_Shape[idx_mag:(idx_mag+ntheta)]

    nfreeparams = num_amplitudes + nspecz
    params_initial = np.zeros(nfreeparams) 
    
    # Optimisation is VERY sensitive to initial values.
    # Set initial alphas to no-mag case (this initialisation seems to achieve best chi^2).
    params_initial[num_amplitudes:] = np.ones(nspecz) 
else:
    # If not accounting for magnification, just fitting for the amplitude parameter(s)
    nfreeparams = num_amplitudes
    params_initial = np.zeros(nfreeparams)



ndof = (ntomo * nspecz * ntheta) - nfreeparams
result = minimize(chi2, params_initial, args=(gt, cov_inv), options={'maxiter': 200000}, method='BFGS')


Use_Curve_Fit = True # Try using curve_fit instead of chi^2 optimize
                     # This is ONLY for Single_Bin=True
def model_curve_fit(x, amplitude):
    return Dls_over_Ds * amplitude / x

if Use_Curve_Fit and Single_Bin:
    params_0_cf, params_0_cferr = curve_fit(model_curve_fit, thetas, gt, p0=params_initial, sigma=cov )  
    # Save the outputs
    np.savetxt(OUTDIR+'/%sx%s_FitParamsCF_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
                                                                Include_Magnification,Include_IA,
                                                                save_tag,nofz_shift), params_0_cf)
    np.savetxt(OUTDIR+'/%sx%s_FitParamsErrCF_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
                                                                   Include_Magnification,Include_IA,
                                                                   save_tag,nofz_shift), params_0_cferr)
    np.savetxt(OUTDIR+'/%sx%s_FitModelCF_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
                                                               Include_Magnification,Include_IA,
                                                               save_tag,nofz_shift),
               np.transpose(np.vstack(( thetas,func(params_0_cf) )) ) )


    
chi2_red = result['fun']/ndof
p_value = 1 - stats.chi2.cdf(result['fun'], ndof)

print("chi^2=%1.2f, dof=%i, chi^2/dof=%1.2f, p=%1.3f, success=%s\n" % (result['fun'], ndof, chi2_red, p_value, result['success']))

f = open(OUTDIR+'/%sx%s_Shear_ratio%s.res'%(SOURCE_TYPE,LENS_TYPE,save_tag), 'w')
f.write("chi^2=%1.2f, dof=%i, chi^2/dof=%1.2f, p=%1.3f, success=%s\n" % (result['fun'], ndof, chi2_red, p_value, result['success']))
f.close()

params_0=result['x']
# Based on the following stackoverflow page, it looks like
# the sqrt of the inverse Hessian matrix diag is the error on the fitted parameters
# https://stackoverflow.com/questions/43593592/errors-to-fit-parameters-of-scipy-optimize
params_0err=np.sqrt( np.diag(result['hess_inv']) )
#np.savetxt(OUTDIR+'/%sx%s_FitParams_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
#                                                          Include_Magnification,Include_IA,
#                                                          save_tag,nofz_shift), params_0)
#np.savetxt(OUTDIR+'/%sx%s_FitParamsErr_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
#                                                             Include_Magnification,Include_IA,
#                                                             save_tag,nofz_shift), params_0err)
#np.savetxt(OUTDIR+'/%sx%s_FitModel_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
#                                                         Include_Magnification,Include_IA,
#                                                         save_tag,nofz_shift),
#           np.transpose(np.vstack(( thetas,func(params_0) )) ) )
# Save the p-value to be read in and plotted by the code Compare_BFParams_And_Models.py
#np.savetxt(OUTDIR+'/%sx%s_pvalue_Mag%s_IA%s%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,
#                                                       Include_Magnification,Include_IA,
#                                                 save_tag,nofz_shift), np.c_[p_value])

######### plots ###

def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)

fig, axs = plt.subplots(nrows=ntomo, ncols=nspecz, sharex=True, sharey=True)
###adjustFigAspect(fig,aspect=.5)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)

for i in range(ntomo):
    for j in range(nspecz):

        # For plotting we separate the data back into the different source-lens bin combinations
        minplt=i*nspecz*ntheta+ j*ntheta
        maxplt=minplt + ntheta

        thetas_to_plot = thetas[minplt:maxplt]
        gt_to_plot = gt[minplt:maxplt]
        gterr_to_plot = diag[minplt:maxplt]
        #gterr_to_plot = gterr[minplt:maxplt]
        #to do the Dls/Ds functions!!!
        func_to_plot = func(params_0)[minplt:maxplt]
        if Single_Bin:
            # In this case, axs is not subscriptable (it's 1x1 in size)
            tmp_axs = axs      # avoids annoying plotting error
        else:
            tmp_axs = axs[i,j]
        tmp_axs.set_xscale('log')
        tmp_axs.set_xlim(theta_min*0.99,theta_max*2.)
        tmp_axs.set_ylim(-0.0075,0.025)

        tmp_axs.errorbar(thetas_to_plot, thetas_to_plot*gt_to_plot, yerr=thetas_to_plot*gterr_to_plot, fmt='o', capsize=2, markersize=4.)
        tmp_axs.plot(thetas_to_plot, thetas_to_plot*func_to_plot)

        tmp_axs.text(theta_min*2., 0.013, "t "+str(tomo_bins[i]+1)+", sp "+str(spec_bins[j]+1))
        tmp_axs.axhline(0., color='gray', linestyle=':', linewidth=1.)

if Single_Bin:
    axs.set_ylabel(r"$\theta \times <\gamma_t>$ [arcmin]")
    axs.set_xlabel(r"$\theta$ [arcmin]")
else:
    axs[2,0].set_ylabel(r"$\theta \times <\gamma_t>$ [arcmin]")
    for i in range(nspecz):
        axs[-1,i].set_xlabel(r"$\theta$ [arcmin]")
        
fig.suptitle(str(LENS_TYPE)+ " lenses, " + str(SOURCE_TYPE) +" sources,  $\chi^2/$dof={:1.2f}".format(chi2_red)+",  $p$={:2.1f}%".format(p_value*100.) )
plt.savefig(OUTDIR+'/%sx%s_Shear_ratio%s%s%s.png'%(SOURCE_TYPE,LENS_TYPE,DFLAG,save_tag,nofz_shift))
#plt.show()

