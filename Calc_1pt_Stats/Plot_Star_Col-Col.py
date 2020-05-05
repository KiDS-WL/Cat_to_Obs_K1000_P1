# 12/02/2020, B. Giblin, Postdoc, Edinburgh
# Produce a colour-colour diagram for stars in KiDS


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
from scipy.stats import binned_statistic_2d
import sys
import time
# The following are used to cross-match (RA,Dec)'s between
# the sources and PSF data catalogues:
from astropy.coordinates import SkyCoord
from astropy import units as u

#from fitting import * # MV scripts                                                                                                
# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}
plt.rc('font', **font)


Read_Cat_Or_Pickle = "Pickle" # Set to "Cat" to read data from individual fields
                           # (~890 s for the 492 North fields)
                           # Set to "Pickle" to read in pre-pickled file
                           # (~15s for the 492 North fields)


# EDIT: 05/05/2020
# The SG_FLAG=0 objects in the source catalogue contains ~7.5% galaxy contamination even after colour cuts.
# (These contaminants are removed in Lensfit prior to PSF modelling).
# But to remove them here, we need to cross-match the objects in the source catalogue
# with those in the PSF data catalouge FOR EACH FIELD INDIVIDUALLY, BEFORE APPLYING colour-colour and SG_FLAG cuts.
LFver = '321'
expname = 'PSFRES_XI_glab_%s' %LFver # Needed to locate relevant PSF exposure catalogs

                           
# Define the Baldry locus
def Baldry(gmi):
    if (gmi < 0.3):
        flocus = -0.7172
    elif ((gmi<2.3) & (gmi>=0.3)) :
        flocus = -0.89 + 0.615*gmi - 0.13*gmi**2
    else :
        flocus = -0.1632
    return flocus + 0.07  # need to add 0.07 to the Baldry criteria to match the stellar locus                     



# Read in names of KiDS fields
def Get_Field_Names(NorS):
    fname = '/home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/K1000_%s.txt' %NorS
    with open(fname) as f:
        Fields_tmp = f.readlines()

    Fields=[]
    for x in Fields_tmp:
        Fields.append(x.split()[1])

    return Fields

# CODE UP: Read in South and stack the two?
Fields = Get_Field_Names("N")
Cycle_Array = Fields[0:]      # How many fields to use?

if Read_Cat_Or_Pickle == "Cat": # READ IN DATA BY CYCLING THROUGH INDIVIDUAL FIELDS

    # Cycle through fields reading in relevant data
    data_indir = "/disk09/KIDS/KIDSCOLLAB_V1.0.0/" # Contains all the KiDS Field directories
    # ARRAYS TO STORE THE DATA
    MAG_GAAP_g = np.zeros([])
    MAG_GAAP_i = np.zeros([])
    MAG_GAAP_J = np.zeros([])
    MAG_GAAP_Ks = np.zeros([])
    SG_FLAG = np.zeros([])
    RA = np.zeros([])
    Dec= np.zeros([])

    t1 = time.time()
    i = 0
    for f in Cycle_Array:
        print("Reading in field %s of %s"%(i,len(Cycle_Array)))

        # ----- Loop through the 5 PSF catalogues corresponding to the 5 exposures ----- #
        psf_exposures = glob.glob('%s/%s/checkplots/%s/*_PSFres.cat'%(data_indir,f,expname))
        tmp_PSF_RA = np. zeros([])
        tmp_PSF_Dec = np.zeros([])
        for e in psf_exposures:
            psf_fitsfile = fits.open(e)
            tmp_PSF_RA = np.append(tmp_PSF_RA, psf_fitsfile[1].data['ALPHA_J2000'] )
            tmp_PSF_Dec =np.append(tmp_PSF_Dec,psf_fitsfile[1].data['DELTA_J2000'] )
        # Delete the 1st elements which are from initiliastion. 
        tmp_PSF_RA = np.delete(tmp_PSF_RA, 0)
        tmp_PSF_Dec = np.delete(tmp_PSF_Dec, 0)
        # ------------------------------------------------------------------------------- #
        
            
        # Now read the corresponding data from the source catalogues
        # Think there should only be one of these files per Field
        # (not one per exposure). But just in case, loop:
        exposures = glob.glob('%s/%s/r_SDSS/colourcat_V1.0.0A/KIDS_*_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask.cat'%(data_indir,f))
        for e in exposures:
            fitsfile = fits.open(e)
            tmp_RA = fitsfile[1].data['RAJ2000']
            tmp_Dec= fitsfile[1].data['DECJ2000']
            # Cross-match the tmp_(RA,Dec)'s with the tmp_PSF_(RA,Dec)'s
            c = SkyCoord(ra=tmp_PSF_RA*u.degree, dec=tmp_PSF_Dec*u.degree)
            catalog = SkyCoord(ra=tmp_RA*u.degree, dec=tmp_Dec*u.degree)
            idx_cross,_,_ = c.match_to_catalog_sky(catalog) # Returns the elems of tmp_(RA,Dec)
                                                            # which are closest to those in tmp_PSF_(RA,Dec)
            idx_cross = np.unique(idx_cross)                # Remove duplicates: this cross-matching cuts tmp_(RA,Dec) down by factors>100
            
            # Selection criteria to apply
            flag_m = fitsfile[1].data['MASK'][idx_cross]
            flag_g = fitsfile[1].data['FLAG_GAAP_g'][idx_cross]
            flag_i = fitsfile[1].data['FLAG_GAAP_i'][idx_cross]
            flag_J = fitsfile[1].data['FLAG_GAAP_J'][idx_cross]
            flag_Ks = fitsfile[1].data['FLAG_GAAP_Ks'][idx_cross]
            idx_colour = ((flag_m==0) & (flag_g==0) & (flag_i==0) & (flag_J==0) & (flag_Ks==0))

            # colours
            tmp_MAG_GAAP_g = fitsfile[1].data['MAG_GAAP_g'][idx_cross][idx_colour]
            tmp_MAG_GAAP_i = fitsfile[1].data['MAG_GAAP_i'][idx_cross][idx_colour]
            tmp_MAG_GAAP_J = fitsfile[1].data['MAG_GAAP_J'][idx_cross][idx_colour]
            tmp_MAG_GAAP_Ks = fitsfile[1].data['MAG_GAAP_Ks'][idx_cross][idx_colour]
            tmp_SG_FLAG = fitsfile[1].data['SG_FLAG'][idx_cross][idx_colour]
            tmp_MAG_AUTO = fitsfile[1].data['MAG_AUTO'][idx_cross][idx_colour]
            # Apply the cuts to the (RA,Dec) too:
            tmp_RA = tmp_RA[idx_cross][idx_colour]
            tmp_Dec= tmp_Dec[idx_cross][idx_colour]

            # Further get rid of the mags with unphysical values of \pm 99
            # These are non-detections
            idx_unphys = ((abs(tmp_MAG_GAAP_g)<99.) & (abs(tmp_MAG_GAAP_i)<99.) &
                    (abs(tmp_MAG_GAAP_J)<99.) & (abs(tmp_MAG_GAAP_Ks)<99.) &
                    (tmp_MAG_AUTO>18) & (tmp_MAG_AUTO<22.5)) 
        
            MAG_GAAP_g = np.append(MAG_GAAP_g, tmp_MAG_GAAP_g[idx_unphys])
            MAG_GAAP_i = np.append(MAG_GAAP_i, tmp_MAG_GAAP_i[idx_unphys])
            MAG_GAAP_J = np.append(MAG_GAAP_J, tmp_MAG_GAAP_J[idx_unphys])
            MAG_GAAP_Ks = np.append(MAG_GAAP_Ks, tmp_MAG_GAAP_Ks[idx_unphys])
            SG_FLAG = np.append(SG_FLAG, tmp_SG_FLAG[idx_unphys])
            RA = np.append(RA, tmp_RA[idx_unphys])
            Dec = np.append(Dec, tmp_Dec[idx_unphys])
        i+=1
    
    MAG_GAAP_g = np.delete(MAG_GAAP_g, 0)
    MAG_GAAP_i = np.delete(MAG_GAAP_i, 0)
    MAG_GAAP_J = np.delete(MAG_GAAP_J, 0)
    MAG_GAAP_Ks = np.delete(MAG_GAAP_Ks, 0)
    SG_FLAG = np.delete(SG_FLAG, 0)
    RA = np.delete(RA, 0)
    Dec = np.delete(Dec, 0)

    t2 = time.time()
    print( "It took %.f s to read in the colour info for %s Fields" %((t2-t1),len(Cycle_Array)) )
    
    # Pickle Data
    print("Pickling the data from %s fields" %len(Cycle_Array) )
    t1 = time.time()
    np.save( 'MAG_giJKs_SGFLAG_%sFields_CrossMatched'%len(Cycle_Array),
            np.column_stack(( MAG_GAAP_g, MAG_GAAP_i, MAG_GAAP_J, MAG_GAAP_Ks, SG_FLAG, RA, Dec )) )
    print("Pickling the data from %s Fields took %.1f s." %( len(Cycle_Array),(t2-t1)) )

elif Read_Cat_Or_Pickle == "Pickle":
    
    print( "Reading in pre-pickled data for %s Fields." %len(Cycle_Array) )
    t1 = time.time()
    data = np.load('MAG_giJKs_SGFLAG_%sFields_CrossMatched.npy' %len(Cycle_Array))
    MAG_GAAP_g = data[:,0]
    MAG_GAAP_i = data[:,1]
    MAG_GAAP_J = data[:,2]
    MAG_GAAP_Ks = data[:,3]
    SG_FLAG = data[:,4]
    
    t2 = time.time()
    print( "Reading in pre-pickled data for %s Fields took %.1f s." %(len(Cycle_Array),(t2-t1)) )

    

# Function to make a density plot in the
# (MAG_GAAP_g - MAG_GAAP_i) VS (MAG_GAAP_J - MAG_GAAP_Ks) parameter space
def MeanQ_VS_XY(Q, X,Y, num_XY_bins): 
    sumQ_grid, yedges, xedges, binnum = binned_statistic_2d(Y, X, Q, statistic='sum', bins=num_XY_bins)
    AvQ_grid = sumQ_grid #/ np.sum(sumQ_grid) # normalise
    return AvQ_grid, yedges, xedges


# Mag-band differences
gmi = (MAG_GAAP_g - MAG_GAAP_i)
JmK = (MAG_GAAP_J - MAG_GAAP_Ks)
# Some of these have extreme values - narrow the range
idx_gmi = np.where(np.logical_and(gmi>-1, gmi<4))[0]
idx_JmK = np.where(np.logical_and(JmK>-1.5, JmK<1.5))[0]
idx_cut = np.intersect1d( idx_gmi, idx_JmK ) # Get the common elements
gmi_cut = gmi[idx_cut]
JmK_cut = JmK[idx_cut]
# Now make the density plot in colour-space
nbins = 100    # dimensions of output grid
cgrid, JmK_bins, gmi_bins = MeanQ_VS_XY(np.ones_like(gmi_cut),
                                        gmi_cut, JmK_cut, nbins)
                                        

# Also create a density plot for the SG_FLAG=0 objects
SG_FLAG_cut = SG_FLAG[ idx_cut ]
gmi_cut_SG0 = gmi_cut[ SG_FLAG_cut == 0 ]   # Of those in narrowed colour range...
JmK_cut_SG0 = JmK_cut[ SG_FLAG_cut == 0 ]   # ...find those with SG_FLAG=0
SG0_grid, SG0_JmK_bins, SG0_gmi_bins = MeanQ_VS_XY(np.ones_like(gmi_cut_SG0),
                                                   gmi_cut_SG0, JmK_cut_SG0, nbins)

def Find_Height_At_Contour_Area(Grid, Area):
    sorted_L = np.sort(np.ravel(Grid))          # flatten to 1D array and put in order
    sorted_L = sorted_L[::-1]                   # put max elements first-
                                                # want to start at peak of likelihood.                                                                                 
    Total_Prob = np.sum(sorted_L)
    # scroll down from the peak of the distribution: adding the heights = adding the volume                                 # since each column has the same base area. Stop when you've reached queried Area
    for i in range(len(sorted_L)):
        if np.sum(sorted_L[:i]) / Total_Prob > Area:
            Height = sorted_L[i-1]
            break
    return Height
cgrid_contours = [ Find_Height_At_Contour_Area(cgrid, 0.99),
                   Find_Height_At_Contour_Area(cgrid, 0.95),
                   Find_Height_At_Contour_Area(cgrid, 0.68) ]
SG0_grid_contours = [ Find_Height_At_Contour_Area(SG0_grid, 0.99),
                      Find_Height_At_Contour_Area(SG0_grid, 0.95),
                      Find_Height_At_Contour_Area(SG0_grid, 0.68) ]


# Plot the colour-colour density plot
def Plot_CC_Grid(grid1,xvals1,yvals1,levels1,
                 grid2,xvals2,yvals2,levels2,
                 scatter_x, scatter_y,
                 title):
    
    #cmap = plt.cm.binary
    cmap = plt.cm.rainbow
    #cmap.set_bad('w', 1.)
    fig = plt.figure(figsize=(10.5,7))
    ax = fig.add_subplot(111)
        
    im = ax.imshow(grid1, cmap=cmap, interpolation='none',
                    vmin=grid1.min(), vmax=grid1.max(), origin='lower', 
                    extent=[np.min(xvals1), np.max(xvals1), np.min(yvals1), np.max(yvals1)],
                    aspect='equal')
    ##plt.contour(xvals1, yvals1, grid1, levels=levels1, colors='darkblue')
    
    # Scatter the points outside the outer contour
    ax.scatter( scatter_x, scatter_y, color='magenta', s=1, marker='.')
    ax.contour(xvals2, yvals2, grid2, levels=levels2, colors='magenta', linewidths=3)

    ax.plot(gmi_loc,JmK_loc, color='r', linewidth=3)
    ax.plot(gmi_loc,JmK_loc+0.2,color='r', linewidth=3)

    ax.set_xlim([xvals1.min(), xvals1.max()])
    ax.set_ylim([yvals1.min(), yvals1.max()])
    ax.set_xlabel(r'$g-i$') 
    ax.set_ylabel(r'$J-K_{\rm{s}}$') 
    #ax.set_title(title)
    ax.set_aspect('auto')

    # Shift subplot to add colorbar on the right-hand side
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical', label=r'Number of sources $[\times 10^3]$')
    #plt.colorbar()
    plt.savefig('star-col-col_%sFields.png' %len(Cycle_Array))
    plt.show() 
    return 

# Extract the SG0 objects outside of the outer contour
# so you can separately scatter them

# First convert the SG0 gmi & JmK arrays into bins on the grid:
bins_x = ( nbins * (gmi_cut_SG0 - gmi_cut_SG0.min()) / (gmi_cut_SG0.max() - gmi_cut_SG0.min()) ).astype(int)
bins_x[np.argmax(bins_x)] += -1    # Correct the one (max) bin on edge of grid range
bins_y = ( nbins * (JmK_cut_SG0 - JmK_cut_SG0.min()) / (JmK_cut_SG0.max() - JmK_cut_SG0.min()) ).astype(int)
bins_y[np.argmax(bins_y)] += -1
# Return the idx of those outside of the outer contour:
idx_outside = np.where( SG0_grid[bins_y, bins_x] < SG0_grid_contours[0] )[0]
scatter_x_outside = gmi_cut_SG0[ idx_outside ]
scatter_y_outside = JmK_cut_SG0[ idx_outside ]


#overplot Baldry locus
gmi_loc = np.arange(-1.0,4.0,0.1)
# define vectorized Baldry
Baldry_v = np.vectorize(Baldry)
JmK_loc = Baldry_v(gmi_loc)

# now determine the fraction of galaxies in the star sample with delta_sg_JK > 0.2                                 
delta_sg_JK = JmK_cut_SG0 - Baldry_v(gmi_cut_SG0)
N_tot = float( len(delta_sg_JK) )
N_contam = len(delta_sg_JK[delta_sg_JK>0.2])

print ("Fraction of galaxies in star sample:", N_contam/N_tot)

Plot_CC_Grid(cgrid/1e3, gmi_bins[:-1],JmK_bins[:-1],cgrid_contours,
             SG0_grid, SG0_gmi_bins[:-1], SG0_JmK_bins[:-1],SG0_grid_contours,
             scatter_x_outside, scatter_y_outside,
             '%s Fields, SG_FLAG=0'%len(Cycle_Array))
