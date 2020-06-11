import numpy as np
import astropy.io.fits as fits

nBins  = 5
nTheta = 9 
nPairs = int(nBins*(nBins+1)/2)
blind='C'


filename="../../../kids/fits/xipm_KIDS1000_Blind"+blind+"_no_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits"

F=fits.open(filename)
ext=F["COVMAT"]
covariance= ext.data

filename='../../../kids/multiplicative_bias/m_cov_r_1.00.ascii'
file  = open(filename)
m_cov = np.loadtxt(file)

cov_m_all = np.zeros((2*nTheta*nPairs,2*nTheta*nPairs))

# where plus ends and minus starts
plus = nPairs * nTheta

for z1 in range(nBins):
	for z2 in range(z1, nBins):
		for z3 in range(z1, nBins):
			if(z1==z3):
				start_z4 = z2
			else:
				start_z4 = z3
			for z4 in range(start_z4,nBins):
				Cm1= m_cov[z1,z3]
				Cm2= m_cov[z2,z4]
				Cm3= m_cov[z1,z4]
				Cm4= m_cov[z2,z3]
				Cm_tot=(Cm1+Cm2+Cm3+Cm4)
				p1=int(nBins*(z1)-(z1)*(z1-1)/2+(z2)-(z1))
				p2=int(nBins*(z3)-(z3)*(z3-1)/2+(z4)-(z3))
				# print(z1,z2,' ',z3,z4)
				# print(p1,p2)
# 
				# read in the xipm input values
				filename='../input_for_xipm_sigma_m_covariance/chain/output_test_A/shear_xi_plus_binned/bin_'+str(z2+1)+'_'+str(z1+1)+'.txt'
				file  = open(filename)
				xip_z1_z2 = np.loadtxt(file)
# 
				filename='../input_for_xipm_sigma_m_covariance/chain/output_test_A/shear_xi_plus_binned/bin_'+str(z4+1)+'_'+str(z3+1)+'.txt'
				file  = open(filename)
				xip_z3_z4 = np.loadtxt(file)
# 
				filename='../input_for_xipm_sigma_m_covariance/chain/output_test_A/shear_xi_minus_binned/bin_'+str(z2+1)+'_'+str(z1+1)+'.txt'
				file  = open(filename)
				xim_z1_z2 = np.loadtxt(file)
# 
				filename='../input_for_xipm_sigma_m_covariance/chain/output_test_A/shear_xi_minus_binned/bin_'+str(z4+1)+'_'+str(z3+1)+'.txt'
				file  = open(filename)
				xim_z3_z4 = np.loadtxt(file)
#
				Cov_pp = np.outer( xip_z1_z2 , xip_z3_z4 )
				Cov_mm = np.outer( xim_z1_z2 , xim_z3_z4 )
				Cov_pm = np.outer( xip_z1_z2 , xim_z3_z4 )
				Cov_mp = np.outer( xim_z1_z2 , xip_z3_z4 )
				Cov_sig_m_pp = Cm_tot * Cov_pp
				Cov_sig_m_mm = Cm_tot * Cov_mm
				Cov_sig_m_pm = Cm_tot * Cov_pm
				Cov_sig_m_mp = Cm_tot * Cov_mp
# 
				imin=nTheta*p1
				imax=imin+nTheta
				jmin=nTheta*p2
				jmax=jmin+nTheta
				cov_m_all[imin:imax,jmin:jmax] = Cov_sig_m_pp
				cov_m_all[jmin:jmax,imin:imax] = Cov_sig_m_pp
				cov_m_all[imin+plus:imax+plus,jmin:jmax] = Cov_sig_m_mp
				cov_m_all[jmin:jmax,imin+plus:imax+plus] = Cov_sig_m_mp
				cov_m_all[imin:imax,jmin+plus:jmax+plus] = Cov_sig_m_pm
				cov_m_all[jmin+plus:jmax+plus,imin:imax] = Cov_sig_m_pm
				cov_m_all[imin+plus:imax+plus,jmin+plus:jmax+plus] = Cov_sig_m_mm
				cov_m_all[jmin+plus:jmax+plus,imin+plus:imax+plus] = Cov_sig_m_mm

cov_tot = covariance + cov_m_all
filename = '../blind'+blind+'/thps_cov_kids1000_xipm_matrix_with_sigma_m.dat'
np.savetxt(filename,cov_tot)


import matplotlib.pyplot as plt
plt.clf()
fig, ax = plt.subplots()
im = ax.imshow(cov_m_all)
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()

plt.clf()
fig, ax = plt.subplots()
im = ax.imshow(cov_tot)
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()


corr=np.zeros((len(cov_tot),len(cov_tot)))
for i in range(len(cov_tot)):
    for j in range(len(cov_tot)):
        corr[i,j]=cov_tot[i,j]/np.sqrt(cov_tot[i,i]*cov_tot[j,j])

plt.clf()
fig, ax = plt.subplots()
im = ax.imshow(corr)
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()


# filename='../../kids1000_chains/covariance/outputs/Covariance_blindA_nMaximum_20_0.50_300.00_nBins5_sigma_m_from_m_cov_0.0000.ascii'

# file=open(filename)
# cosebis_sigma_m = np.loadtxt(filename)
# plt.clf()
# fig, ax = plt.subplots()
# im = ax.imshow(cosebis_sigma_m)
# cbar = ax.figure.colorbar(im, ax=ax)
# plt.show()

