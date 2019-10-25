#Script from Hendrik - hacked by Catherine
#Now we use treecorr and Mandelbaum estimator
#to do the random and boost correction
#we can simplify this script
#I also removed the scale selection part
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

plt.rc('font', size=10)

ntheta = 4
theta_min = 2.
theta_max = 30.

measurements = {}

thetas_list = []
gt_list = []
gx_list = []
gterr_list = []

tomo_list = []
speczbin_list = []
Dls_over_Ds_list = []

SPINDIR='/disk09/KIDS/K1000_TWO_PT_STATS//OUTSTATS/SHEAR_RATIO/'
ntomo=5
nspecz=5
nspin=1000

for tomobin in range(ntomo):
    for speczbin in range(nspecz):

        gtfile=SPINDIR+'/GT/K1000_GT_6Z_source_'+str(tomobin+1)+'_5Z_lens_'+str(speczbin+1)+'.asc'
        gtdat=ascii.read(gtfile)

        measurements['thetas_'+str(speczbin)+"_"+str(tomobin)] = gtdat['meanr']
        measurements['gt_'+str(speczbin)+"_"+str(tomobin)]     = gtdat['gamT']
        measurements['gx_'+str(speczbin)+"_"+str(tomobin)]     = gtdat['gamX']
        measurements['gterr_'+str(speczbin)+"_"+str(tomobin)]  = gtdat['sigma']

        thetas_list.append(measurements['thetas_'+str(speczbin)+"_"+str(tomobin)])
        gt_list.append(measurements['gt_'+str(speczbin)+"_"+str(tomobin)])
        gx_list.append(measurements['gx_'+str(speczbin)+"_"+str(tomobin)])
        gterr_list.append(measurements['gterr_'+str(speczbin)+"_"+str(tomobin)])

        tomo_list.append(tomobin*np.ones((ntheta),dtype=np.int16))
        speczbin_list.append(speczbin*np.ones((ntheta),dtype=np.int16))

        # Read in the Dls_over_Ds data created with Dls_over_Ds.py
        Dls_over_Ds_file = 'Dls_over_Ds_data/Dls_over_Ds_DIR_6Z_source_'+str(tomobin+1)+'_5Z_lens_'+str(speczbin+1)+'.asc'        
        Dls_over_Ds_tmp = np.loadtxt(Dls_over_Ds_file)
        measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)] = np.repeat(Dls_over_Ds_tmp,ntheta)
        Dls_over_Ds_list.append(measurements['Dls_over_Ds_'+str(speczbin)+"_"+str(tomobin)])
        
thetas = np.hstack(thetas_list)
tomo = np.hstack(tomo_list)
speczbin = np.hstack(speczbin_list)
gt = np.hstack(gt_list) 
gx = np.hstack(gx_list) 
gterr = np.hstack(gterr_list)
Dls_over_Ds = np.hstack(Dls_over_Ds_list)

nmatrix = nspecz*ntomo*ntheta
gtlens=np.zeros([ntheta,nspin])
diag=np.zeros(nmatrix)
covdiag=np.zeros([nmatrix,nmatrix])

for tomobin in range(ntomo):
    for speczbin in range(nspecz):
        for ispin in range(nspin):
            gtfile=SPINDIR+'/GT/SPIN/K1000_GT_SPIN_'+str(ispin)+'_6Z_source_'+str(tomobin+1)+'_5Z_lens_'+str(speczbin+1)+'.asc'
            gtspindat=ascii.read(gtfile)
            gtlens[:,ispin]=gtspindat['gamT'] 
        if speczbin==0 and tomobin==0:      
            gtprev = np.copy(gtlens)
        else:    
            gtall = np.vstack((gtprev,gtlens))
            gtprev= np.copy(gtall)
            print len(gtall)
    
cov=np.cov(gtall)

print len(cov)
for i in range(nmatrix):
    diag[i]=np.sqrt(cov[i,i])
    covdiag[i,i]=cov[i,i]

Corr=np.corrcoef(gtall)
plt.imshow(Corr, interpolation='None')
plt.colorbar()
plt.axis('off')
#plt.show()
plt.savefig('K1000xBOSS_Shear_ratio_correlation_matrix.png')

cov_inv=np.linalg.inv(cov)
#cov_inv=np.linalg.inv(covdiag)

#############################

def func(params):
    amplitude = np.hstack((params,)*ntomo)
    return Dls_over_Ds * amplitude

######## including covariance ###

def chi2(params, data, cov_inv):
    return np.dot(data-func(params),np.dot(cov_inv,data-func(params)))

result = minimize(chi2, np.zeros([ntheta*nspecz]), args=(gt, cov_inv), options={'maxiter': 150000})

ndof = (ntomo * nspecz * ntheta) - (ntheta * nspecz)
chi2_red = result['fun']/ndof
p_value = 1 - stats.chi2.cdf(result['fun'], ndof)

print "chi^2=%1.2f, dof=%i, chi^2/dof=%1.2f, p=%1.3f, success=%s\n" % (result['fun'], ndof, chi2_red, p_value, result['success'])

f = open('K1000xBOSS_Shear_ratio.res', 'w')
f.write("chi^2=%1.2f, dof=%i, chi^2/dof=%1.2f, p=%1.3f, success=%s\n" % (result['fun'], ndof, chi2_red, p_value, result['success']))
f.close()

params_0=result['x']

#np.savetxt(model_file,np.transpose(np.vstack((thetas,func(params_0)))))

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

        axs[i,j].set_xscale('log')
        axs[i,j].set_xlim(theta_min*0.99,theta_max*2.)
        axs[i,j].set_ylim(-0.0075,0.025)

        axs[i,j].errorbar(thetas_to_plot, thetas_to_plot*gt_to_plot, yerr=thetas_to_plot*gterr_to_plot, fmt='o', capsize=2, markersize=4.)
        axs[i,j].plot(thetas_to_plot, thetas_to_plot*func_to_plot)

        axs[i,j].text(theta_min*2., 0.013, "t "+str(i+1)+", sp "+str(j+1))
        axs[i,j].axhline(0., color='gray', linestyle=':', linewidth=1.)


axs[2,0].set_ylabel(r"$\theta \times <\gamma_t>$ [arcmin]")
for i in range(nspecz):
    axs[-1,i].set_xlabel(r"$\theta$ [arcmin]")
fig.suptitle(" BOSS lenses, K1000-P1 sources,  $\chi^2/$dof={:1.2f}".format(chi2_red)+",  $p$={:2.1f}%".format(p_value*100.))
plt.savefig('K1000xBOSS_Shear_ratio.png')
plt.show()
