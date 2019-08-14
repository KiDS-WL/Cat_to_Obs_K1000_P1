#
import numpy as np
import ldac
import matplotlib.pyplot as plt

import sys
from scipy.stats import binned_statistic_2d

from matplotlib import rcParams

#rcParams['ps.useafm'] = True
#rcParams['pdf.use14corefonts'] = True
#rcParams.update({'figure.autolayout': True})
#rcParams['axes.linewidth'] = 2 # set the value globally

#plt.rc('font', family='serif')

#font = {'family' : 'serif',
#        'weight' : 'normal',
#        'size'   : 19}

#plt.rc('font', **font)
#plt.rc('xtick', labelsize=26) 
#plt.rc('ytick', labelsize=26) 

patch=sys.argv[1]
LFVER=sys.argv[2]
infile='%s_PSF_residuals_V1.0.0_%s.cat' %(patch, LFVER)
#infile='%s_PSF_residuals_V1.0.0_svn_309b.cat' %(patch)

# Define custom colormaps: Set pixels with no sources to white
cmap = plt.cm.jet
cmap.set_bad('w', 1.)

cmap_multicolor = plt.cm.jet
cmap_multicolor.set_bad('w', 1.)

# Create figure and subplots
fig = plt.figure()
fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95,
                    bottom=0.07, top=0.95)

#open the ldac catalogue using functions in ldac.py

ldac_cat = ldac.LDACCat(infile)
ldac_table = ldac_cat['OBJECTS']

# read in the useful columns

KiDS_RA=ldac_table['ALPHA_J2000']
KiDS_Dec=ldac_table['DELTA_J2000']
#Xpos=ldac_table['Xpos_THELI']
#Ypos=ldac_table['Ypos_THELI']
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

print size_data
print size_mod

print min(chip), max(chip)

for j in range (1,2):   # use 6,7 for star count

    if(j==1):

        de1 = pe1_data - pe1_mod
        de2 = pe2_data - pe2_mod
        
    elif(j==2):

        de1 = pe1_data
        de2 = pe2_data  

    elif(j==3):
        de1 = pe1_mod
        de2 = pe2_mod

    elif(j==4):

        de1 = size_data - size_mod
        de2 = size_data - size_mod
        
    else:

        de1 = size_data
        de2 = size_mod          

    de1_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de1,statistic='mean', bins=20, range=[[0.0,21000.],[0.0,21000.0]])
    de2_mean, yedges, xedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic='mean', bins=20, range=[[0.0,21000.],[0.0,21000.0]])
    star_count, ysedges, xsedges, binnum = binned_statistic_2d(Ypos, Xpos, de2,statistic='count', bins=50, range=[[0.0,21000.],[0.0,21000.0]])

#    de1_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de1,statistic='mean', bins=[32,10])
#    de2_mean, xedges, yedges, binnum = binned_statistic_2d(chip, mag, de2,statistic='mean', bins=[32,10])

    for i in range(1,3):  # use 3,4 for star count, 1,3 for other plots

        if(i==1):
            plt.subplot(223)
            bindat = de1_mean

            print bindat
        elif(i==2):
            plt.subplot(224)
            bindat = de2_mean

        elif(i==3):
            bindat = star_count  # if you want to look at the number density
            xedges=xsedges
            yedges=ysedges

            print xedges, yedges
            
        im = plt.imshow(bindat, origin='lower',
                        extent=[xedges[0], xedges[20],
                                yedges[0], yedges[20]],
                        aspect='equal', interpolation='none', cmap=cmap)
        #im = plt.imshow(bindat, origin='lower',
        #                extent=[yedges[0], yedges[10],
        #                        xedges[0], xedges[32]],
        #                aspect='auto', interpolation='none', cmap=cmap)

        #plt.xlim(xedges[0], xedges[10])
        #plt.ylim(yedges[0], yedges[31])
        plt.gca().set_aspect('equal',adjustable='box')
        plt.xlabel('Xpos')
        plt.ylabel('Ypos')

        #plt.xlabel('star mag')
        #plt.ylabel('chip_ID')

        if(j==1):
            if(i==1):
                plt.clim(-0.002,0.002)
                plt.title('depsf 1')
            else:
                plt.clim(-0.002,0.002)
                plt.title('depsf 2')
        elif(j==2):
            plt.clim(-0.015,0.015)
            if(i==1):
                plt.title('PSF e1 (data)')
            else:
                plt.title('PSF e2 (data)')
        elif(j==3):
            plt.clim(-0.015,0.015)
            if(i==1):
                plt.title('PSF e1 (model)')
            else:
                plt.title('PSF e2 (model)')
        elif(j==4):
            plt.clim(-0.03,0.03)
            if(i==1):
                plt.title('dsize (data - model)')
            else:
                plt.title('dsize (data - model)')
        elif(j==5):
            plt.clim(1.9,2.4)
            if(i==1):
                plt.title('PSF size (data)')
            else:
                plt.title('PSF size (model)')      
        else:
            plt.clim(0.0,1000)  # star count
            plt.title('Nstars')                
                    
        plt.gca().set_aspect('equal', adjustable='box')
        #plt.gca().set_aspect('auto',adjustable='box',anchor='N')
        plt.tight_layout()
        
    # Make an axis for the colorbar on the right side

    cax = fig.add_axes([0.1, 0.05, 0.8, 0.01])
    fig.colorbar(im, cax=cax, orientation='horizontal')

    if(j==1):
        plt.savefig('PLOTS/depsf_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))
    elif(j==2):
        plt.savefig('PLOTS/epsf_data_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))
    elif(j==3):
        plt.savefig('PLOTS/epsf_model_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))
    elif(j==4):
        plt.savefig('PLOTS/desize_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))
    elif(j==5):
        plt.savefig('PLOTS/size_data_model_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))
    else:
        plt.savefig('PLOTS/star_count_%s_V1.0.0_%s_X_Y.png' % (patch, LFVER))


    plt.show()
