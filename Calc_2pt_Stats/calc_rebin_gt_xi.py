# ----------------------------------------------------------------
# File Name:           calc_rebin_gt_xi.py
# Author:              Chieh-An Lin (calin@roe.ac.uk)
# Description:         
# ----------------------------------------------------------------

import treecorr
import sys
import numpy as np


def rebin(r_min, r_max, N_r, lin_not_log, meanr, meanlnr, weight, valueBlock, wgtBlock):
    if lin_not_log == 'true':
        bAndC  = np.linspace(r_min, r_max, 2*N_r+1)
    else:
        bAndC  = np.logspace(np.log10(r_min), np.log10(r_max), 2*N_r+1)
    
    ctrBin    = bAndC[1::2] ## [arcmin]
    bins      = bAndC[0::2]
    wgt_r     = weight * meanr
    wgt_lnr   = weight * meanlnr
    wgt_value = weight * valueBlock
    N_col_v   = valueBlock.shape[0]
    N_col_w   = wgtBlock.shape[0]
    
    if wgt_value.ndim > 1:
        wgt_value = wgt_value.T
    if wgtBlock.ndim > 1:
        wgtBlock  = wgtBlock.T
    
    
    binned_r          = []
    binned_lnr        = []
    binned_valueBlock = []
    binned_wgtBlock   = []
    
    for i in range(N_r):
        ind = (meanr > bins[i]) * (meanr < bins[i+1])
        
        if ind.any():
            wgt_sum = weight[ind].sum()
            binned_r.append(wgt_r[ind].sum() / wgt_sum)
            binned_lnr.append(wgt_lnr[ind].sum() / wgt_sum)
            binned_valueBlock.append(wgt_value[ind].sum(axis=0) / wgt_sum)
            binned_wgtBlock.append(wgtBlock[ind].sum(axis=0))
        else:
            print("WARNING: not enough bins to rebin to "+str(N_r)+" log bins")
            binned_r.append(np.nan)
            binned_lnr.append(np.nan)
            binned_valueBlock.append([np.nan]*N_col_v)
            binned_wgtBlock.append([np.nan]*N_col_w)
    
    binned_r   = np.array(binned_r)
    binned_lnr = np.array(binned_lnr)
    binned_valueBlock = np.array(binned_valueBlock).T
    binned_wgtBlock   = np.array(binned_wgtBlock).T
    
    return ctrBin, binned_r, binned_lnr, binned_valueBlock, binned_wgtBlock


if __name__ == '__main__':
    
    # Read in user input to set
    if len(sys.argv) <7: 
        print("Usage: %s loadName GT_or_XI N_theta theta_min(arcmin) theta_max(arcmin) lin_not_log? saveName" % sys.argv[0]) 
        sys.exit(1)
    else:
        loadName = sys.argv[1]
        GT_or_XI = sys.argv[2]
        N_theta = int(sys.argv[3]) 
        theta_min = float(sys.argv[4]) 
        theta_max = float(sys.argv[5]) 
        lin_not_log = sys.argv[6]
        saveName = sys.argv[7]

    ## Read old TreeCorr file
    data = np.loadtxt(loadName, comments='#', delimiter=None).T

    if GT_or_XI == 'XI':
        meanr = data[1]
        meanlnr = data[2]
        weight = data[9]

        valueBlock = data[3:9]
        wgtBlock = data[9:11]

        ## Turn sigma into sigma^2
        valueBlock[4:] = valueBlock[4:]**2

        ## Rebin
        ctrBin, binned_r, binned_lnr, binned_valueBlock, binned_wgtBlock = rebin(theta_min, theta_max, N_theta, lin_not_log, meanr, meanlnr, weight, valueBlock, wgtBlock)

        ## Turn sigma^2 into sigma
        valueBlock[4:] = np.sqrt(valueBlock[4:])

        # Write it out to a file and praise-be for Jarvis and his well documented code
        treecorr.util.gen_write(saveName,
            ['r_nom', 'meanr', 'meanlogr', 'xip', 'xim', 'xip_im', 'xim_im', 'sigma_xip', 'sigma_xim', 'weight', 'npairs'],
            [ctrBin, binned_r, binned_lnr, 
            binned_valueBlock[0], binned_valueBlock[1], binned_valueBlock[2], binned_valueBlock[3], binned_valueBlock[4], binned_valueBlock[5],
            binned_wgtBlock[0], binned_wgtBlock[1]],
            precision=12
        )
    
    else:
        meanr = data[1]
        meanlnr = data[2]
        weight = data[6]

        indList = [3, 4, 5, 8, 9, 10, 11, 12]
        valueBlock = data[indList]
        
        indList = [6, 7]
        wgtBlock = data[indList]

        ## Turn sigma into sigma^2
        indList = [2, 7]
        valueBlock[indList] = valueBlock[indList]**2

        ## Rebin
        ctrBin, binned_r, binned_lnr, binned_valueBlock, binned_wgtBlock = rebin(theta_min, theta_max, N_theta, lin_not_log, meanr, meanlnr, weight, valueBlock, wgtBlock)

        ## Turn sigma^2 into sigma
        valueBlock[indList] = np.sqrt(valueBlock[indList])

        # Write it out to a file and praise-be for Jarvis and his well documented code
        treecorr.util.gen_write(saveName,
            ['r_nom','meanr','meanlogr','gamT','gamX','sigma','weight','npairs', 'nocor_gamT', 'nocor_gamX', 'rangamT','rangamX','ransigma' ],
            [ctrBin, binned_r, binned_lnr, 
              binned_valueBlock[0], binned_valueBlock[1], binned_valueBlock[2],
              binned_wgtBlock[0], binned_wgtBlock[1], 
              binned_valueBlock[3], binned_valueBlock[4], binned_valueBlock[5], binned_valueBlock[6], binned_valueBlock[7]],
            precision=12
        )


