import sys
import numpy as np
import cosmology
import astropy.io.fits as fits
#From Hendrik
# 25/05/2020: Edited by B. Giblin,
# now includes the delta-z shift and error to the nofz
# if a third input argument containing "_nofzUpp" or "_nofzLow" is input.

#from astropy.cosmology import WMAP9 as cosmo
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from scipy import integrate

LCDM=cosmology.Cosmology()

lens_cat = fits.open(sys.argv[1])
lens_data = lens_cat[1].data
z_lens = lens_data['Z']

z_lens_max = 2.0
no_z_lens_bins = 200
z_lens_step = z_lens_max / no_z_lens_bins

lens_hist=np.histogram(z_lens,range=(0.,z_lens_max),bins=no_z_lens_bins,normed=True)

norm_lens = np.sum(lens_hist[0])

lens_hist_norm = lens_hist[0] / norm_lens

# Annoyingly, the K1000 source data binned using SOM-z's are FITS files,
# whereas the DIR-binned K1000 data & MICE data are ascii catalogues,
# necessitating the following if statement:
if '.fits' in sys.argv[2]:
    source_cat = sys.argv[2]
    tmp = fits.open(source_cat)
    zs = tmp[1].data['binstart']
    nzs = tmp[1].data['density']
    source_hist = np.column_stack(( zs, nzs ))
else:
    source_hist = np.loadtxt(sys.argv[2])

# Apply the delta-z shifts and errors if
# a third input argument says to.
try:
    nofz_shift=sys.argv[3]
except IndexError:
    nofz_shift=""

if "_nofz" in nofz_shift:
    # Read in delta-z shifts and errors
    dz,dz_err = np.loadtxt('nofz_shifts_perbin_mean_err.asc', usecols=(0,1), unpack=True)
    tbin = int(source_cat.split('TOMO')[-1].split('_')[0]) # Extracted tomo bin number

    # Here we figure out what kind of shift is applied:
    # Up means we find the n(z-A*dz_err), and shift it UP to redshift z
    # Down means we find n(z+A*dz_err), and shift DOWN to redshift z.
    # Mix means we shift EVEN tomo bins UP and ODD tomo bins DOWN.

    shift_size = int( nofz_shift.split('sig')[0][-1] ) # How many sigma shifts to apply to n(z)
    if "Up" in nofz_shift:
        shift_size *= -1.

    # If "Down" is in nofz_shift, then we leave shift_size as it is (+ve number).
    
    elif "Mix" in nofz_shift:
        # Then we're shifting UP if it's an even bin, or DOWN if it's an odd bin:
        if (tbin-1) % 2 ==0:
            shift_size *= -1.

    # Apply the shift:
    new_zs = source_hist[:,0] + shift_size * dz_err[tbin-1] # forget the mean shift # +shift_size*dz[tbin-1]
    if "Up" in nofz_shift:
        # Get rid of the ~1 element that goes negative from this shift
        new_zs[new_zs<0.] = 0.

    # Interp old nofz to new shifted redshifts
    new_nzs = np.interp( new_zs, source_hist[:,0], source_hist[:,1] )
    # Replace old source histogram with the shifted version
    source_hist[:,1] = new_nzs
        
z_source_max = source_hist[-1,0]

no_z_source_bins = np.shape(source_hist)[0]
z_source_step = source_hist[1,0] - source_hist[0,0]

norm_source = np.sum(source_hist[:,1])

source_hist_norm = source_hist[:,1] / norm_source

def Dls_over_Ds(z1,z2):
    if z2>z1:
        return LCDM.dang(z1,z2)/LCDM.dang(0.,z2)
    else:
        return 0.

integral=0.

for i in range(no_z_lens_bins):
    for j in range(no_z_source_bins):
        if (i+0.5)*z_lens_step < (j+0.5)*z_source_step:
            integral += lens_hist_norm[i] * source_hist_norm[j] * Dls_over_Ds((i+0.5)*z_lens_step,(j+0.5)*z_source_step)

print(integral)

#### Checking astropy for consistency (slower) ###
#
#def Dls_over_Ds2(z1,z2):
#    return cosmo.angular_diameter_distance_z1z2(z1,z2)/cosmo.angular_diameter_distance(z2)
#
#integral=0.
#
#for i in range(no_z_lens_bins):
#    for j in range(no_z_source_bins):
#        if (i+0.5)*z_lens_step < (j+0.5)*z_source_step:
#            integral += lens_hist_norm[i] * source_hist_norm[j] * Dls_over_Ds2((i+0.5)*z_lens_step,(j+0.5)*z_source_step)
#
#print(integral)

#### Checking against direct numerical integration (not working yet) ###
#
#def dn_dz_lens(z):
#    if z <= z_lens_max:
#        z_index = int(z / z_lens_step)
#        return lens_hist_norm[z_index] * no_z_lens_bins / z_lens_max
#    else:
#        return 0.
#
#def dn_dz_source(z):
#    if z < z_source_max:
#        z_index = int(z / z_source_step)
#        return source_hist_norm[z_index] * no_z_source_bins / z_source_max
#    else:
#        return 0.
#
#def func(z1,z2):
#    return dn_dz_lens(z1) * dn_dz_source(z2) * Dls_over_Ds(z1,z2)
#
#print(integrate.dblquad(func, 0., 10., lambda x: 0., lambda x: 10.))
