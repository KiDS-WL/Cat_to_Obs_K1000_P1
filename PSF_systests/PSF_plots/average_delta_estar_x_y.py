#
import numpy as np
import ldac
import matplotlib.pyplot as plt
import sys
from scipy.stats import binned_statistic_2d
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rcParams

#============================================================

# Input information
patch=sys.argv[1]   # e.g K1000_N or K1000_S or K1000
LFVER=sys.argv[2]   # e.g glab_321
XY_or_chip=sys.argv[3] # X_Y or chip

# Input residual file created with paste_psf_residual_cats.sh
infile='%s_PSF_residuals_V1.0.0_%s.cat' %(patch, LFVER)

#============================================================
# Set up the figure
fig = plt.figure(figsize=(8, 6))

# Define custom colormaps: Set pixels with no sources to white
cmap = plt.cm.jet
cmap.set_bad('w', 1.)

# Set the image grid
if(XY_or_chip=='X_Y'):  # square
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(2,2),
                 axes_pad=0.25,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="edge",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
elif(XY_or_chip=='chip'): # row
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,4),
                 axes_pad=0.25,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="edge",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )
else:
        print "Bin in X/Y or chip/mag?  Define XY_or_chip as X_Y or chip"
        exit

# Some font setting
rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}

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

for j in range (1,3):  # what do we want to bin?  And we may want to fix the colour scale limits (cmin,cmax)

    if(j==1):   # data

        de1 = pe1_data 
        de2 = pe2_data
        cmin=-0.015
        cmax= 0.015
        labinfo = '$\epsilon^*$'

    elif(j==2):  # ellipticity residuals

        de1 = pe1_data - pe1_mod
        de2 = pe2_data - pe2_mod 
        cmin=-0.001
        cmax= 0.001
        labinfo = '$\delta\epsilon^*$'

#    elif(j==3): # you might also want to look at size residuals and size
                 # but to make pretty publication plots I comment this part out

#        de1 = size_data - size_mod
#        de2 = size_mod   

    if(XY_or_chip=='X_Y'):

        de1_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de1,statistic='mean', bins=20, range=[[0.0,21000.],[0.0,21000.0]])
        de2_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic='mean', bins=20, range=[[0.0,21000.],[0.0,21000.0]])
        star_count, ysedges, xsedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic='count', bins=20, range=[[0.0,21000.],[0.0,21000.0]])

    elif(XY_or_chip=='chip'):

        de1_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de1,statistic='mean', bins=[32,10])
        de2_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de2,statistic='mean', bins=[32,10])
        star_count, xsedges, ysedges, binnum = binned_statistic_2d(chip, mag, de2,statistic='count', bins=[32,10])


    for i in range(1,3):  # now make the figures
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
            im = ax.imshow(bindat, origin='lower',
                        extent=[xedges[0], xedges[20],
                                yedges[0], yedges[20]],
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
grid[0].set_title('$\epsilon_1^*$',size=15)
grid[1].set_title('$\epsilon_2^*$',size=15)

if(XY_or_chip=='X_Y'):
    grid[0].set_ylabel('Y',size=14)     
    grid[2].set_ylabel('Y',size=14)   
    grid[2].set_xlabel('X',size=14)     
    grid[3].set_xlabel('X',size=14) 
else:
    grid[0].set_ylabel('Chip ID',size=14) 
    grid[2].set_xlabel('$r$-band magnitude            ',size=14) 
    grid[2].set_title('$\delta\epsilon_1^*$',size=15)
    grid[3].set_title('$\delta\epsilon_2^*$',size=15)

plt.savefig('PLOTS/epsf_and_depsf_%s_V1.0.0_%s_%s.png' % (patch, LFVER, XY_or_chip))
plt.show()
