#
import numpy as np
import ldac
import matplotlib.pyplot as plt
#import pylab as plt
import sys
from scipy.stats import binned_statistic_2d
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rcParams
from matplotlib import rc

#============================================================

# Input information
patch=sys.argv[1]   # e.g K1000_N or K1000_S or K1000
LFVER=sys.argv[2]   # e.g glab_321
XY_or_chip=sys.argv[3] # X_Y or chip

# Input residual file created with paste_psf_residual_cats.sh
indir = '/home/cech/KiDSLenS/Cat_to_Obs_K1000_P1/PSF_systests/PSF_plots/' # Line added by Giblin to access CH's PSF catalogue.
infile='%s/%s_PSF_residuals_V1.0.0_%s.cat' %(indir, patch, LFVER)

#============================================================
# Define custom colormaps: Set pixels with no sources to white
cmap = plt.cm.jet
cmap.set_bad('w', 1.)


# Set the figure and image grid
if(XY_or_chip=='X_Y'):  # square
    fig = plt.figure(figsize=(6, 6))    # was (8,6)
    nrows=3                             # Changed to 3 to include sigma_ePSF
    ncols=2
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(nrows,ncols),
                 axes_pad=0.25,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="edge",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
elif(XY_or_chip=='chip'): # row
    fig = plt.figure(figsize=(6, 6))
    nrows=1
    ncols=4
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(nrows,ncols),
                 axes_pad=0.25,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="edge",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
else:
        print ( "Bin in X/Y or chip/mag?  Define XY_or_chip as X_Y or chip" )
        exit

# Some font setting
# NOTE: Where the following line appears (here, or after setting of font size) has funny effects on fontsize.
rcParams["mathtext.fontset"] = "cm" # This makes episilons appear as epsilons rather than varepsilons 
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

#plt.rcParams["mathtext.fontset"] = "cm" # This makes episilons appear as epsilons rather than varepsilons  
plt.rc('font', **font)

#============================================================
#Now read in the data
#open the ldac catalogue using functions in ldac.py

ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the useful columns



KiDS_RA=ldac_table['ALPHA_J2000']
KiDS_Dec=ldac_table['DELTA_J2000']
Xpos=ldac_table['Xpos']
Ypos=ldac_table['Ypos']

chip=ldac_table['chip_ID']
mag=ldac_table['STAR_MAG']

pe1_data=ldac_table['psfe_data0']
pe2_data=ldac_table['psfe_data1']

pe1_mod=ldac_table['psfe_model0']
pe2_mod=ldac_table['psfe_model1']

Q11=ldac_table['moments_data0']
Q22=ldac_table['moments_data1']
Q12=ldac_table['moments_data2']

size_data= np.sqrt(Q11*Q22-Q12**2)

Q11=ldac_table['moments_model0']
Q22=ldac_table['moments_model1']
Q12=ldac_table['moments_model2']

size_mod= np.sqrt(Q11*Q22-Q12**2)

#============================================================
#Now make different sub plots
# count the grid square that I want to plot in
gridpos=-1

for j in range (1,nrows+1):  # what do we want to bin?  And we may want to fix the colour scale limits (cmin,cmax)

    if(j==1):   # <ePSF>

        de1 = pe1_data 
        de2 = pe2_data
        cmin=-0.015
        cmax= 0.015
        labinfo = r'$\overline{ \epsilon^{\rm PSF} }$'
        stat='mean'

    elif(j==2): # sigma[ ePSF ]
        de1 = pe1_data
        de2 = pe2_data
        cmin= 0.014
        cmax= 0.030
        labinfo = r'$\sigma_{\epsilon^{\rm PSF}}$'
        stat='std'
        
    elif(j==3):  # <delta ePSF>

        de1 = pe1_data - pe1_mod
        de2 = pe2_data - pe2_mod 
        cmin=-0.001
        cmax= 0.001
        labinfo = r'$\overline{ \delta\epsilon^{\rm PSF} }$'
        stat='mean'
        
#    elif(j==3): # you might also want to look at size residuals and size
                 # but to make pretty publication plots I comment this part out

#        de1 = size_data - size_mod
#        de2 = size_mod   

    if(XY_or_chip=='X_Y'):

        de1_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de1,statistic=stat, bins=20, range=[[0.0,21000.],[0.0,21000.0]])
        de2_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic=stat, bins=20, range=[[0.0,21000.],[0.0,21000.0]])
        star_count, ysedges, xsedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic='count', bins=20, range=[[0.0,21000.],[0.0,21000.0]])

        # In std case, get zeros instead on nans at edge:
        if stat == 'std':
            de1_mean[ de1_mean==0. ] = np.nan
            de2_mean[ de2_mean==0. ] = np.nan
        

    elif(XY_or_chip=='chip'):

        de1_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de1,statistic=stat, bins=[32,10], range=[[0.5,32.5],[18.0,22.6]])
        de2_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de2,statistic=stat, bins=[32,10], range=[[0.5,32.5],[18.0,22.6]])
        star_count, xsedges, ysedges, binnum = binned_statistic_2d(chip, mag, de2,statistic='count', bins=[32,10])


    for i in range(1,ncols+1):  # now make the figures
        gridpos=gridpos + 1

        if(i==1):
            bindat = de1_mean

        elif(i==2):
            bindat = de2_mean



        #elif(i==3):
        #    bindat = star_count  # if you want to look at the number density
                                  # to make pretty publication plots I comment this part out  
        #    xedges=xsedges
        #    yedges=ysedges

        ax=grid[gridpos]

        if(XY_or_chip=='X_Y'):
            scale_pix=1000.
            im = ax.imshow(bindat, origin='lower',
                        extent=[xedges[0]/scale_pix, xedges[20]/scale_pix,
                                yedges[0]/scale_pix, yedges[20]/scale_pix],
                        aspect='equal', interpolation='none', cmap=cmap, clim=(cmin,cmax))


        elif(XY_or_chip=='chip'):             
            im = ax.imshow(bindat, origin='lower',
                        extent=[yedges[0], yedges[10],
                                xedges[0], xedges[32]],
                        aspect='equal',interpolation='none', cmap=cmap, clim=(cmin,cmax)) 

        thecb = ax.cax.colorbar(im)  
        thecb.set_label_text(labinfo,size=15)

        # lets make the axes a little bit thicker
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.25)
        ax.tick_params(width=1.25)

# And fix the labels with a nice large font size 
grid[0].set_title(r'$\epsilon_1^{\rm PSF}$',size=15)
grid[1].set_title(r'$\epsilon_2^{\rm PSF}$',size=15)

if(XY_or_chip=='X_Y'):
    grid[0].set_ylabel(r'$Y \, [\times 10^3 \, \rm{pixels}]$',size=14)     
    grid[2].set_ylabel(r'$Y \, [\times 10^3 \, \rm{pixels}]$',size=14)
    grid[4].set_ylabel(r'$Y \, [\times 10^3 \, \rm{pixels}]$',size=14)
    
    grid[4].set_xlabel(r'$X \, [\times 10^3 \, \rm{pixels}]$',size=14)     
    grid[5].set_xlabel(r'$X \, [\times 10^3 \, \rm{pixels}]$',size=14) 
else:
    grid[0].set_ylabel('Chip ID',size=14) 
    grid[2].set_xlabel('$r$-band magnitude            ',size=14) 
    grid[2].set_title('$\delta\epsilon_1^*$',size=15)
    grid[3].set_title('$\delta\epsilon_2^*$',size=15)

plt.savefig('PLOTS/epsf_and_depsf_%s_V1.0.0_%s_%s.png' % (patch, LFVER, XY_or_chip))
plt.show()
