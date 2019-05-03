# 01/05/2019, B. M. Giblin, PhD Student, Edinburgh
# Code that accepts a FITS file containing e_1,2 & m_1,2, and calculates optimal weights
# and then saves another FITS file with this info in it.

import numpy as np
import sys
from astropy.io import fits
import pylab as plt

# READ USER INPUT FILENAME
userinput = sys.argv[1:]
filename = userinput[0]
fitsfile = fits.open(filename)

# Stuff that's currently soft-coded here, but could be given as input?
nbins = 8													# Average quantities in this many bins and interpolate between them
savename = filename.split('.fits')[0] + '_OptW.fits'		# You can change this savename to whatever
															# Currently it just pegs '_OptW.fits to the end of the input filename.
selection = fitsfile[1].data['fitclass'] == 0.				# Object selection criteria


# Pull the quantities with column names 'x_string' and 'y_string' from the fitsfile
# and apply object selection specified by 'selection'
def Pull_FITS(fitsfile, x_string, y_string, selection):
	yvalues = (fitsfile[1].data[ y_string ])[ selection ]	
	xvalues = (fitsfile[1].data[x_string])[ selection ]
	return xvalues, yvalues



# Take either weighted-average OR the sum of the ydata within the nbins bins of xdata
def Binning(nbins, xdata, ydata, weights, AvgOrSum):

	x_dist, bins_x_edges = np.histogram(xdata, bins=nbins)
	width = bins_x_edges[1] - bins_x_edges[0]
	bins_x= bins_x_edges[:-1] + width/2.
	y_dist = np.zeros(nbins)
	y_dist_err = np.zeros_like(y_dist)
	for i in range(nbins):
  	  idx = np.where( np.logical_and(xdata>bins_x_edges[i], xdata<bins_x_edges[i+1]) )[0]
  	  if len(idx) != 0:
			y_dist[i] = np.sum( weights[idx]*ydata[idx] ) 
			y_dist_err[i] = np.std( ydata[idx] ) 		# GIBLIN THIS ERROR IS NOT RIGHT CURRENTLY
														# REMIND YERSEN HOW TO DO WEIGHTED STAN-DEV.

			if AvgOrSum == "Avg":
				y_dist[i] /= np.sum( weights[idx] )
				y_dist_err[i] /= np.sum( weights[idx] ) 

	return bins_x, y_dist, y_dist_err



def Plot_xVSy(x, y, y_err, titlestring):
	plt.figure()
	plt.errorbar( x, y, yerr=y_err, fmt='o', markersize=8, color='magenta' )
	plt.title(titlestring)
	plt.show()
	return



def Calc_OptWeights(fitsfile, nbins, selection):

	# KEY:
	# '_pg' means 'per galaxy, '_pmb' means 'per magnitude bin'


	# First calculate the unweighted variance
	mag_pg, e1_pg = Pull_FITS(fitsfile, 'mag_out', 'e1', selection)
	mag_pg, e2_pg = Pull_FITS(fitsfile, 'mag_out', 'e2', selection)
	mag_pg, e1c_pg = Pull_FITS(fitsfile, 'mag_out', 'e1_corr', selection)
	mag_pg, e2c_pg = Pull_FITS(fitsfile, 'mag_out', 'e2_corr', selection)
	e_sq_pg = (e1_pg-e1c_pg)**2 + (e2_pg-e2c_pg)**2
	mag_pmb, var_e_pmb, var_e_err_pmb = Binning(nbins, mag_pg, e_sq_pg, np.ones_like(e_sq_pg), "Avg")
	var_e_pmb *= 0.5 			# Factor 1/2 is not included in the weighted avg'ing done by Binning function
	var_e_err_pmb *= 0.5
	Plot_xVSy(mag_pmb, np.sqrt(var_e_pmb), var_e_err_pmb/(2*np.sqrt(var_e_pmb)), r'$\sigma_e$ VS mag')

	# The optimal weight per magnitude bin
	mag_pg, m1_pg = Pull_FITS(fitsfile, 'mag_out', 'metacal_m1', selection)		# We checked, makes no diff. if you use m1 or m2.
	dummy, m1_pmb, m1_err_pmb = Binning(nbins, mag_pg, m1_pg, np.ones_like(m1_pg), "Avg")

	weight_pmb = (1+m1_pmb) / var_e_pmb 
	weight_err_pmb = weight_pmb * np.sqrt( (m1_err_pmb/m1_pmb)**2 + (var_e_err_pmb/var_e_pmb)**2 )
	Plot_xVSy(mag_pmb, weight_pmb, weight_err_pmb, 'Optimal weight VS mag')

	# Interpolate this to get an optimal weight per galaxy
	weight_pg = np.interp( mag_pg, mag_pmb, weight_pmb)
	weight_err_pg = np.interp( mag_pg, mag_pmb, weight_err_pmb)

	return mag_pmb, weight_pmb, weight_err_pmb, weight_pg, weight_err_pg

mag_pmb, weight_pmb, weight_err_pmb, weight_pg, weight_err_pg = Calc_OptWeights(fitsfile, nbins, selection)




# Re-save the FITS file under 'newfitsname' with an extra column of data in it.
def Save_NewCol_2Fits(fitsfile, data, colname, newfitsname):
	print("Saving new FITS file with new column-name %s under %s" %(colname,newfitsname))
	print("This may take a wee minute...")
	cols = [] 
	cols.append(fits.Column(name=colname, format='D', array=data))
	orig_cols = fitsfile[1].columns
	new_cols = fits.ColDefs(cols)
	new_fitsfile = fits.BinTableHDU.from_columns(orig_cols + new_cols)
	try:
		new_fitsfile.writeto(newfitsname)
	except IOError:
		print("There already is a file saved as %s" %newfitsname)
		usr_input = raw_input("Do you want to overwrite? Type y for yes, anything else will be treated as a no...")
		if usr_input == "y":
			print("Okay, overwriting...")
			new_fitsfile.writeto(newfitsname, clobber=True)
		else:
			print("Okay, not overwriting.")		
	return


# Save a new fitsfile with the optimal weights as an extra column

# First need to get the optimal for EVERY galaxy, not just the 'selection' ones
selection_all = []
for i in range(len(selection)):
	selection_all.append(True)
mag_pg_all, dummy = Pull_FITS(fitsfile, 'mag_out', 'e1', selection_all)
weight_pg_all = np.interp( mag_pg_all, mag_pmb, weight_pmb)

Save_NewCol_2Fits(fitsfile, weight_pg_all, 'Optimal_LFweight', savename)


	





