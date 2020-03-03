# ----------------------------------------------------------------
# File Name:           calc_gt_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate gamma_t
#                      we're using the Mandelbaum estimator where the randoms are subtracted
#                      so we need a lens, source and random catalogue 
#                      script will need to change if keywords in KIDS cats are updated
#  Treecorr doc: https://rmjarvis.github.io/TreeCorr/_build/html/correlation2.html
# ----------------------------------------------------------------

import treecorr
import sys
import numpy as np
import astropy.io.fits as fits

## Calculate shape noise given noisy & noise-free shear
def subtractNoise(g_1, g_2, eps_1, eps_2):
    g   = g_1 + 1j * g_2
    g_c = g_1 - 1j * g_2
    eps = eps_1 + 1j * eps_2
    
    e = (eps - g) / (1.0 - g_c*eps)
    e = np.array([e.real, e.imag])
    return e

if __name__ == '__main__':
    # Read in user input to set the nbins, theta_min, theta_max, lin_not_log, lenscat, rancat, sourcecat, outfilename, weighted
    if len(sys.argv) < 10: 
        print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) lin_not_log(true or false)? lenscat.fits \
               randomcat.fits sourcecat.fits outfilename weighted_analysis?" % sys.argv[0]) 
        sys.exit(1)
    else:
        nbins = int(sys.argv[1]) 
        theta_min = float(sys.argv[2]) 
        theta_max = float(sys.argv[3]) 
        lin_not_log = sys.argv[4] 
        lenscatname = sys.argv[5]
        rancatname = sys.argv[6]
        sourcecatname = sys.argv[7]
        outfile = sys.argv[8]
        weighted = sys.argv[9]

    # prepare the catalogues

    ## For mocks
    if 'MFP_galCat' in lenscatname:
        
        ## Lens catalogue
        L1_data = fits.getdata(lenscatname, 1)
        L1_RA   = L1_data.field('ALPHA_J2000')
        L1_DEC  = L1_data.field('DELTA_J2000')
        L1_wgt  = L1_data.field('weight')
        
        ## Random catalogue
        R1_data = fits.getdata(rancatname, 1)
        R1_RA   = R1_data.field('ALPHA_J2000')
        R1_DEC  = R1_data.field('DELTA_J2000')
        
        ## Get the weights for randoms
        if '2dFLenS' in rancatname or 'MFP_selection' in rancatname:
            R1_wgt  = R1_data.field('WEICOMP')
        else:
            R1_wgt  = R1_data.field('WEIMAG') ## Dummy

        ## Mask selection
        if '2dFLenS' in rancatname or 'BOSS' in rancatname:
            KIDSMASK = R1_data.field('KIDSMASK')
            bitmask  = 0x4000 ## wcs-like
            ind      = ~(KIDSMASK == bitmask)
            R1_RA    = R1_RA[ind]
            R1_DEC   = R1_DEC[ind]
            R1_wgt   = R1_wgt[ind]

        ## Source catalogue
        S2_data = fits.getdata(sourcecatname, 1)
        S2_RA   = S2_data.field('ALPHA_J2000')
        S2_DEC  = S2_data.field('DELTA_J2000')
        S2_g1   = S2_data.field('g1')
        S2_g2   = S2_data.field('g2')
        S2_eps1 = S2_data.field('e1')
        S2_eps2 = S2_data.field('e2')
        S2_wgt  = S2_data.field('weight')

        ## Make TreeCorr instances for lenses & randoms
        lenscat = treecorr.Catalog(ra=L1_RA, dec=L1_DEC, w=L1_wgt, ra_units='deg', dec_units='deg')
        rancat = treecorr.Catalog(ra=R1_RA, dec=R1_DEC, w=R1_wgt, ra_units='deg', dec_units='deg')
        
        ## Make a TreeCorr instance for sources
        if 'obs' in outfile:
            sourcecat = treecorr.Catalog(ra=S2_RA, dec=S2_DEC, g1=S2_eps1, g2=S2_eps2, w=S2_wgt, ra_units='deg', dec_units='deg')
        elif 'signal' in outfile:
            sourcecat = treecorr.Catalog(ra=S2_RA, dec=S2_DEC, g1=S2_g1, g2=S2_g2, w=S2_wgt, ra_units='deg', dec_units='deg')
        elif 'noShear' in outfile:
            S2_e = subtractNoise(S2_g1, S2_g2, S2_eps1, S2_eps2)
            sourcecat = treecorr.Catalog(ra=S2_RA, dec=S2_DEC, g1=S2_e[0], g2=S2_e[1], w=S2_wgt, ra_units='deg', dec_units='deg')
        else:
            raise ValueError('Houston, we have a problem.')

    ## For data
    ## Things get very complicated that this section needs to be revised
    else:
        lenscat = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', 
                                   ra_units='deg', dec_units='deg', w_col='WEICOMP')
        #TODO Catherine, you might want to check if the column names are still different between
        #TODO BOSS & 2dFLenS
        if '2dFLenS' in rancatname:
            w_col = 'WEICOMP'
        else:
            w_col = 'WEIMAG'
        rancat = treecorr.Catalog(rancatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', 
                                  ra_units='deg', dec_units='deg', w_col=w_col)
        sourcecat = treecorr.Catalog(sourcecatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', 
                                     ra_units='deg', dec_units='deg', g1_col='e1', g2_col='e2', w_col='weight')


    #TODO to do include wcs cut
    #,flag_col='KIDSMASK',ignore_flag=16384

    ## Set bin_slop
    if nbins > 100: ## Fine-binning
        inbinslop_NN = 1.0 # when using fine bins I find this is suitable
        inbinslop_NG = 1.2 
    else: ## Broad bins
        inbinslop_NN = 0.03
        inbinslop_NG = 0.05

    # Define the binning based on command line input
    if(lin_not_log=='true'): 
        # Number of source lens pairs
        nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_type='Linear', bin_slop=inbinslop_NN)
        # Number of source random pairs     
        nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_type='Linear', bin_slop=inbinslop_NN)
        # Average shear around lenses     
        ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_type='Linear', bin_slop=inbinslop_NG)
        # Average shear around randoms     
        rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_type='Linear', bin_slop=inbinslop_NG)

    else:
        # Set-up the different correlations that we want to measure

        # Number of source lens pairs
        nlns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_slop=inbinslop_NN)
        # Number of source random pairs     
        nrns = treecorr.NNCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_slop=inbinslop_NN)
        # Average shear around lenses     
        ls = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_slop=inbinslop_NG)
        # Average shear around randoms     
        rs = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
             bin_slop=inbinslop_NG)

    ## A twig from Linc; he likes to use only 8 processors
    if 'MFP_galCat' in lenscatname:
        num_threads = 8
    else:
        num_threads = None

    # Now calculate the different 2pt correlation functions
    nlns.process(lenscat, sourcecat, num_threads=num_threads)
    nrns.process(rancat, sourcecat, num_threads=num_threads)
    ls.process(lenscat, sourcecat, num_threads=num_threads)
    rs.process(rancat, sourcecat, num_threads=num_threads)

    # We will use the Mandelbaum 2006 estimator which includes both the random and boost correction.
    # It is given by
    # gt = (SD-SR)/NsNr
    # SD = sum of shear around source-lens pairs
    # SR = sum of shear around source-random pairs
    # NrNs = number of source random pairs
    # Note that we have a different number of randoms from lenses so we also need to make
    # sure that we rescale NrNs accordingly

    # The easiest way to do this in Treecorr is
    # gt = (SD/NlNs)*NlNs/NrNs - SR/NrNs
    # where NsNl = number of source lens pairs
    # 
    # SD/NsNl = ls.xi
    # NlNs = nlns.weight/nlns.tot
    # NrNs = nrns.weight/nrns.tot
    # SR/NrNs = rs.xi

    gamma_t = ls.xi*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi
    gamma_x = ls.xi_im*(nlns.weight/nlns.tot)/(nrns.weight/nrns.tot) - rs.xi_im

    if (weighted=='true'):   
        # prepare the weighted_square catalogues - hack so that Treecorr returns the correct Npairs for a weighted sample

        lenscat_wsq = treecorr.Catalog(lenscatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                       w_col='WEICOMPsq')
        sourcecat_wsq = treecorr.Catalog(sourcecatname, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                       g1_col='e1', g2_col='e2', w_col='weightsq')

        # Define the binning based on command line input
        if(lin_not_log=='true'): 
            # Average shear around lenses     
            ls_wsq = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                            bin_type='Linear', bin_slop=inbinslop_NG)
        else: # Log is the default bin_type for Treecorr   
            ls_wsq = treecorr.NGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                            bin_slop=inbinslop_NG)

        # Calculate the weighted square 2pt correlation function
        ls_wsq.process(lenscat_wsq, sourcecat_wsq)
        # Calculate the weighted Npairs = sum(weight_a*weight_b)^2 / sum(weight_a^2*weight_b^2)

        npairs_weighted = (ls.weight)*(ls.weight)/ls_wsq.weight

        #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
        treecorr.util.gen_write(outfile,
                ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs_weighted', 'nocor_gamT', 'nocor_gamX', 
                'rangamT','rangamX','ransigma' ],
                [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), npairs_weighted, ls.npairs,
                ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ],
                precision=12)

    else:     
         #Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
         treecorr.util.gen_write(outfile,
                ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 
                'rangamT','rangamX','ransigma' ],
                [ ls.rnom,ls.meanr, ls.meanlogr,gamma_t, gamma_x, np.sqrt(ls.varxi), ls.weight, ls.npairs,
                ls.xi, ls.xi_im, rs.xi, rs.xi_im, np.sqrt(rs.varxi) ],
                precision=12)

