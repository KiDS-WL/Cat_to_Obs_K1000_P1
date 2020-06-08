# 07/04/2020, B. M. Giblin, Postdoc, Edinburgh
# Read in and compare the model fits:

Compare_Mag = False         # Compare fit params/model including/excluding
                           # the impact of magnification.
Compare_Single_Bin = True  # Compare fit params obtained for each l+s bin individually
                           # across lens bins (should be consistent)

import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.gridspec as gridspec

plt.rc('font', size=18)
# ----- Load Input ----- #                                                                                                               
from Get_Input import Get_Input
paramfile = sys.argv[1]   # e.g. params_KiDSdata.dat                                                                                     
GI = Get_Input(paramfile)
SOURCE_TYPE = GI.Source_Type()
LENS_TYPE = GI.Lens_Type()
RANDOM_TYPE = GI.Random_Type()

Cov_Method = "Spin"   # The method for calculating the gamma_t realisations for use in covariance estimation
INDIR='Output/SOURCE-%s_LENS-%s' %(SOURCE_TYPE, LENS_TYPE)
DlsDIR='Dls_over_Ds_data/SOURCE-%s_LENS-%s' %(SOURCE_TYPE, LENS_TYPE)
if "MICE2" in SOURCE_TYPE:
    # Additional identifiers when working with MICE
    # True/Estimated P(z) and Magnification on/off
    Mag_OnOff = GI.Mag_OnOff()
    Pz = GI.Pz_TrueEstimated()
    SN = GI.SN()
    INDIR += '_Pz%s_SN%s_mag%s' %(Pz,SN,Mag_OnOff)
    DlsDIR += '_Pz%s_SN%s_mag%s' %(Pz,SN,Mag_OnOff)
    
elif "K1000" in SOURCE_TYPE:
    Blind = GI.Blind()
    SOMFLAGNAME = GI.SOMFLAGNAME()
    INDIR += '_Blind%s_SOM%s' %(Blind,SOMFLAGNAME)
    DlsDIR += '_Blind%s_SOM%s' %(Blind,SOMFLAGNAME)
    
OUTDIR = INDIR + '/' + Cov_Method.upper()

nz_source = 5
nz_lens = 5
nz_bins = nz_source*nz_lens
ntheta = 4


if Compare_Mag:
    # Compare the fit params and model for Mag On with Mag Off
    params_m = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True))
    params_nm = np.loadtxt(OUTDIR+'/%sx%s_FitParams_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False))
    thetas, model_m = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,True),
                                 usecols=(0,1), unpack=True)
    thetas, model_nm = np.loadtxt(OUTDIR+'/%sx%s_FitModel_Mag%s.dat'%(SOURCE_TYPE,LENS_TYPE,False),
                                  usecols=(0,1), unpack=True)


    
if Compare_Single_Bin:
    nofz_shift=["_nofzDown", "_nofzUp"]      # Use the SRT results for when the Dls/Ds values were calculated from an nofz
                                             # which has been shifted up ('_nofzUp'), down ('_nofzDown') by (delta_z+delta_z_err)

    # Compare the fit params obtained for each l+s bin, across lens bins...
    params_sl = np.zeros([nz_source, nz_lens])
    params_sl_err = np.zeros_like( params_sl )
    params_sl_shift = np.zeros([len(nofz_shift), nz_source, nz_lens])
    params_sl_shift_err = np.zeros_like( params_sl_shift )

    # ...and compare the best-fit models to the data...
    model_sl = np.zeros([ nz_source, nz_lens, ntheta ])
    model_sl_shift = np.zeros([ len(nofz_shift), nz_source, nz_lens, ntheta ])
    data_sl = np.zeros_like( model_sl )
    # ...including also the p-values and data-point error bars:
    p_values = np.zeros([nz_source, nz_lens])
    p_values_shift = np.zeros([len(nofz_shift), nz_source, nz_lens])
    cov = np.load(OUTDIR+'/GTCovMat_%sx%s_6Z_source_5Z_lens_mCovTrue.npy'%(SOURCE_TYPE,LENS_TYPE))
    
    for i in range(nz_source):
        for j in range(nz_lens):
            data_sl[i,j,:] = np.loadtxt('%s/GT_6Z_source_%s_5Z_lens_%s.asc' %(INDIR,(i+1),(j+1)), usecols=(3,), unpack=True)
            params_sl[i,j] = np.loadtxt(OUTDIR+'/%sx%s_FitParamsCF_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j))
            params_sl_err[i,j] = np.loadtxt(OUTDIR+'/%sx%s_FitParamsErrCF_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j))
            theta, model_sl[i,j,:] = np.loadtxt(OUTDIR+'/%sx%s_FitModelCF_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j),
                                                usecols=(0,1), unpack=True)
            p_values[i,j] = np.loadtxt(OUTDIR+'/%sx%s_pvalue_MagFalse_t%ss%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j))
            # Read in the shifted results:
            for k in range(len(nofz_shift)):
                params_sl_shift[k,i,j] = np.loadtxt(OUTDIR+'/%sx%s_FitParamsCF_MagFalse_t%ss%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j,nofz_shift[k]))
                params_sl_shift_err[k,i,j] = np.loadtxt(OUTDIR+'/%sx%s_FitParamsErrCF_MagFalse_t%ss%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j,nofz_shift[k]))
                model_sl_shift[k,i,j,:] = np.loadtxt(OUTDIR+'/%sx%s_FitModelCF_MagFalse_t%ss%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j,nofz_shift[k]),
                                                usecols=(1,), unpack=True)
                p_values_shift[k,i,j] = np.loadtxt(OUTDIR+'/%sx%s_pvalue_MagFalse_t%ss%s%s.dat'%(SOURCE_TYPE,LENS_TYPE,i,j,nofz_shift[k]))
            
    # Read in the magnification model to overplot
    Magnif_Shape = np.load('/home/bengib/kcap_NewInst/kcap/examples/output_magnification_alpha1.0/SRTparam_Bij.npy').flatten()
    
# ------------------------------------- FUNCTIONS FOR COMPARE MAG ON/OFF ANALYSES -------------------------------------

def Plot_gt_MagOnOff(yaxis_label):
    # Plot the Mag On/Off best-fit models
    fig = plt.figure(figsize = (12,10)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, nz_lens)
    k=0
    for i in range(nz_source):
        for j in range(nz_lens):
            ax1 = plt.subplot(gs[k])
            plt.plot( thetas[k*ntheta:(k+1)*ntheta], thetas[k*ntheta:(k+1)*ntheta]*model_m[k*ntheta:(k+1)*ntheta],
                      color='red', label='Mag On' )
            plt.plot( thetas[k*ntheta:(k+1)*ntheta], thetas[k*ntheta:(k+1)*ntheta]*model_nm[k*ntheta:(k+1)*ntheta],
                      color='blue', label='Mag Off' )

            ax1.set_ylim([-0.005, 0.02])

            if i==nz_source-1:
                # Get an x-axis for bottom-row.
                ax1.set_xlabel(r'$\theta$ [arcmin]')
                if j==nz_lens-1:
                    # Get a legend just for the very bottom-right subplot.
                    ax1.legend(loc='upper left', frameon=False)
            else:
                ax1.set_xticks([])

            if j==0:
                # Get a y-axis for first column
                ax1.set_ylabel(yaxis_label)
            else:
                ax1.set_yticks([])
            k+=1
            ax1.text(0.95,0.95, 't%s, sp%s'%(i+1,j+1), ha="right", va="top", transform=plt.gca().transAxes)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    return

def Plot_MagOnOff_Amplitudes():
    # Compare the amplitudes obtained in the Mag On/Off cases
    
    plt.figure()
    #plt.scatter( np.arange(20)+1, 1e3*params_m[:20], marker='x', color='red', label='Mag On' )
    #plt.scatter( np.arange(20)+1, 1e3*params_nm[:20], marker='o', color='blue', label='Mag Off' )

    plt.scatter( np.arange(20)+1, params_m[:20]/params_nm[:20], marker='o', color='black')
    plt.xlabel( r'bin number' )
    #plt.ylabel( r'Amplitude $[\times 10^{-3}]$' )
    plt.ylabel( r'Amplitude ratio: mag/no-mag' )
    plt.title( r'Amplitude ratio: %.2f pm %.2f' %( np.mean(params_m[:20]/params_nm[:20]),
                                                   np.std( params_m[:20]/params_nm[:20])))
    plt.legend(loc='best', frameon=True)
    plt.show()
    return

if Compare_Mag:
    Plot_gt_MagOnOff(r'$\theta\gamma_{\rm t}^{\rm model}$')
    Plot_MagOnOff_Amplitudes()
    
# ------------------------------------- ^^^FUNCTIONS FOR COMPARE MAG ON/OFF ANALYSES^^^ -------------------------------------





# ------------------------------------- FUNCTIONS FOR COMPARE SINGLE BIN ANALYSIS -------------------------------------

def Model(amplitude, amplitude_err, i, j):
    Dls_over_Ds_file = '%s/Dls_over_Ds_DIR_6Z_source_%s_5Z_lens_%s.asc' %(DlsDIR, i+1, j+1)
    Dls_over_Ds = np.loadtxt(Dls_over_Ds_file)
    model = Dls_over_Ds * amplitude / theta
    model_err = Dls_over_Ds * amplitude_err / theta
    return model, model_err

def Plot_SingleBin_Params():
    MixLog = False # Make some panels with a log scale, some with a lin scale
    log_panels = [2,3,4]
    
    # Define y-limits for each lens bin panel
    lin_ylimits = [ [ 0.0, 0.029], [ 0.01, 0.029], [ -0.029, 0.18], [ -0.1, 0.14], [ -0.45, 0.49] ]
    mix_ylimits = [ [ -0.005, 0.029], [ 0.01, 0.049], [ 0.001, 2.], [ 0.001, 2.], [ 0.001, 2.] ]
    
    # Plot the Mag On/Off best-fit models
    fig = plt.figure(figsize = (8,6)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, 1)
    colors=[ 'red', 'blue', 'orange', 'green', 'magenta' ]
    tomobin_array = np.linspace(0.5, 5.5, 5) #np.arange(nz_source+1)+0.5
    for i in range(nz_source):
        ax1 = plt.subplot(gs[i])
        if MixLog:
            if i in log_panels:
                ax1.set_yscale('log') #, linthreshy=symlogscales[i] )
            ylimits = mix_ylimits[i] 
            save_tag = '-mixlog'
        else:
            ylimits = lin_ylimits[i]
            save_tag = ''

        # Plot the weighted mean and error per lens bin...
        mean = np.average( params_sl[:,i], weights=(1./params_sl_err[:,i]**2) )
        error = np.sqrt( np.sum( params_sl_err[:,i]**2 ) / nz_source )
        #ax1.fill_between( tomobin_array , (mean-error), (mean+error), color='dimgrey', alpha=0.5 )
        #ax1.plot( tomobin_array , (np.zeros_like(tomobin_array)+mean), color='black', linewidth=1 )
        #print( mean, mean-error, mean+error )

        # ------ vvv USE THE SRT RESULTS w/ UP/DOWN SHIFTS IN THE n(z) TO DEFINE AN ERROR ON THE FITTED PARAMS vvv ------

        shift_sys = abs(params_sl_shift[0,:,i] - params_sl_shift[1,:,i] ) # the "systematic" shift in the params due to delta-z, per t-bin
        shift_sys_err = np.sqrt( params_sl_shift_err[0,:,i]**2 + params_sl_shift_err[1,:,i]**2 ) # statistical error on the sys. shift per bin
        overall_shift_err = shift_sys+shift_sys_err
        # Or just try overall MAX range of params:
        max_shift_range = abs( params_sl_shift[0,:,i]+params_sl_shift_err[0,:,i] - (params_sl_shift[1,:,i]-params_sl_shift_err[1,:,i]) )
        shift_mean = np.average( params_sl[2:,i], weights=(1./max_shift_range[2:]**2) )
        shift_mean_err = np.sqrt( np.sum( max_shift_range[2:]**2 ) / 3)  #nz_source)
        ax1.fill_between( tomobin_array , (shift_mean-shift_mean_err), (shift_mean+shift_mean_err), color='yellow', alpha=0.5 )
        ax1.plot( tomobin_array , (np.zeros_like(tomobin_array)+shift_mean), color='black', linewidth=1 )
        #print( error, shift_mean_err )
        
        #ax1.errorbar( np.arange(1,nz_source+1)-0.1, params_sl_shift[0,:,i], yerr=params_sl_shift_err[0,:,i], fmt='v', color=colors[i], edgecolor=None )
        #ax1.errorbar( np.arange(1,nz_source+1)+0.1, params_sl_shift[1,:,i], yerr=params_sl_shift_err[1,:,i], fmt='^', color=colors[i], edgecolor=None )

        # ------ ^^^ USE THE SRT RESULTS w/ UP/DOWN SHIFTS IN THE n(z) TO DEFINE AN ERROR ON THE FITTED PARAMS ^^^ ------ 
        
        # ...plot the smallest error bar for the bin
        #binmin = np.argmin(params_sl_err[:,i])
        #ax1.fill_between( tomobin_array , (params_sl[binmin,i]-params_sl_err[binmin,i]),
        #                  (params_sl[binmin,i]+params_sl_err[binmin,i]), color='dimgrey', alpha=0.5 )
        
        ax1.errorbar( np.arange(1,nz_source+1), params_sl[:,i], yerr=params_sl_err[:,i], fmt='o',
                      color=colors[i], edgecolor=None, solid_capstyle='projecting', capsize=5 )
        ax1.errorbar( np.arange(1,nz_source+1), params_sl[:,i], yerr=max_shift_range, fmt='o',
                      color=colors[i], edgecolor=None, solid_capstyle='projecting', capsize=5 )

        if i==(nz_source-1):
            ax1.set_xlabel(r'Tomographic bin')
            #ax1.set_xticks([])
        else:
            ax1.set_xticks([])

        if i==2:
            ax1.set_ylabel(r'$2\pi (\sigma_{\rm v} /c)^2$')

        ax1.set_ylim([ ylimits[0] , ylimits[1] ])
        #ax1.set_ylim([ 0.5*(params_sl[0,i]-params_sl_err[0,i]).min(),
        #               2.0*(params_sl[0,i]+params_sl_err[0,i]).max() ])
        ax1.text(0.95,0.95, 'sp%s'%(i+1), ha="right", va="top", transform=plt.gca().transAxes)
        
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(OUTDIR+'/%sx%s_FitParams_MagFalse_AllBinsComparison%s.png'%(SOURCE_TYPE,LENS_TYPE,save_tag))
    plt.show()
    return


def Plot_SingleBin_Data_And_Model():
    plt.rc('font', size=16)
    
    # Plot the gamma_t per l&s bin VS the best-fit model when we perform
    # the SRT on each bin individually.
    alpha_boss = 2.0 # ~The level of magnification found for BOSS
    
    fig = plt.figure(figsize = (8,6)) #figsize = (20,14)
    gs = gridspec.GridSpec(nz_source, nz_lens)
    k=0
    for i in range(nz_source):
        for j in range(nz_lens):
            ax1 = plt.subplot(gs[k])
            
            # PLOTTING STUFF
            # Best-fit model
            tmp_model, tmp_model_err = Model(params_sl[i,j], params_sl_err[i,j], i, j)
            plt.errorbar( theta, theta*tmp_model, yerr=tmp_model_err, color='orange', fmt='o')
            plt.plot( theta, theta*model_sl[i,j,:], color='orange', linewidth=2, label=r'Model')

            # Magnification model
            Mag_Model = 2.*(alpha_boss-1.) * Magnif_Shape[k*ntheta:(k+1)*ntheta] #+ model_sl[i,j,:]
            #print( Mag_Model / model_sl[i,j,:] )
            #plt.plot( theta, theta*Mag_Model, color='magenta', linewidth=2, label=r'Magnification')

            # Data
            error = np.sqrt( np.diag( cov[k*ntheta:(k+1)*ntheta, k*ntheta:(k+1)*ntheta] ) )
            plt.errorbar( theta, theta*data_sl[i,j,:], yerr=theta*error, color='blue', fmt='o', label='Data' )

            ax1.set_ylim(-0.0075,0.025)
            ax1.set_xscale('log')
            ax1.set_xlim(1.98, 60.)
            
            if i==nz_source-1:
                # Get an x-axis for bottom-row
                if j==2:
                # In fact, just label the middle column
                    ax1.set_xlabel(r'$\theta$ [arcmin]')
                #if j==nz_lens-1:
                    # Get a legend just for the very bottom-right subplot.
                #    ax1.legend(loc='upper left', frameon=False)
            else:
                ax1.set_xticks([])

            if j>0:
                # kill y-axis for all but the first column
                # Get a y-axis for first column
                ax1.set_yticks([])

            if j==0 and i==2:
                # Only put y-label on first column middle row.
                ax1.set_ylabel(r'$\theta \times \langle \gamma_{\rm t} \rangle \, [ \rm{arcmin} ]$')

            # Certain panels we need to adjust where the text appears or it clashes with the data points
            text_ypos=0.95
            if j==0:
                if i==3 or i==4:
                    text_ypos=0.4
                
            ax1.text(0.95,text_ypos, 't%s, sp%s\n'%(i+1,j+1)+'$p=$%.1f'%(100*p_values[i,j])+r'%',
                     ha="right", va="top", transform=plt.gca().transAxes, fontsize=10)
            k+=1

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(OUTDIR+'/%sx%s_FitModel_MagFalse_AllBinsComparison.png'%(SOURCE_TYPE,LENS_TYPE))
    plt.show()
    return

if Compare_Single_Bin:
    #Plot_SingleBin_Params()
    Plot_SingleBin_Data_And_Model()

# ------------------------------------- ^^^FUNCTIONS FOR COMPARE SINGLE BIN ANALYSIS^^^ -------------------------------------
