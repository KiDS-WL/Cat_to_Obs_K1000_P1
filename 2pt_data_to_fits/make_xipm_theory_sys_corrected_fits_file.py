import astropy.io.fits as fits


FolderNameInputs  = '../data/'
blind = 'A'
cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'
name_tag    = 'with_m_bias' # with_m_bias # no_m_bias # bmodes

filename = FolderNameInputs+'/kids/fits/xipm_sys_corrected_KIDS1000_Blind'+blind+'_'+name_tag+'_'+cat_version+'.fits'
F=fits.open(filename)
ext=F['xip']
xip_data_sys_corrected = ext.data['VALUE']
ext=F['xim']
xim_data_sys_corrected = ext.data['VALUE']
F.close()


filename = FolderNameInputs+'/kids/fits/xipm_KIDS1000_Blind'+blind+'_'+name_tag+'_'+cat_version+'.fits'
F=fits.open(filename)
ext=F['xip']
xip_data = ext.data['VALUE']
ext=F['xim']
xim_data = ext.data['VALUE']
F.close()

xip_sys  = xip_data - xip_data_sys_corrected
xim_sys  = xim_data - xim_data_sys_corrected

filename = FolderNameInputs+'/kids/mock_data/xipm_theory.fits'
F=fits.open(filename)
ext=F['xip']
xip_th = ext.data['VALUE']
ext=F['xim']
xim_th = ext.data['VALUE']
F.close()

xip_th_sys_corrected = xip_th - xip_sys
xim_th_sys_corrected = xim_th - xim_sys

filename = FolderNameInputs+'/kids/mock_data/xipm_theory.fits'
F = fits.open(filename)
ext=F['xip']
data = ext.data
data['VALUE'][:] = xip_th_sys_corrected.copy()
ext=F['xim']
data = ext.data
data['VALUE'][:] = xim_th_sys_corrected.copy()

filename = FolderNameInputs+'/kids/mock_data/xipm_theory_psf_sys_corrected.fits'
F.writeto(filename)

F = fits.open(filename)
ext=F['xip']
xip = ext.data['VALUE']
