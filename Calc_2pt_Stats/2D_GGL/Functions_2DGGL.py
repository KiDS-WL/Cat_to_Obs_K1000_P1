import numpy as np


def Select_Patch(Q, ra, dec, rlo, rhi, dlo, dhi):
    # return the ra-dec data INSIDE the (ra,dec) range
    # Q is an arbitrary quantity with same len as ra,dec
    idx_ra = np.where(np.logical_and(ra<rhi, ra>rlo))[0]
    idx_dec = np.where(np.logical_and(dec<dhi, dec>dlo))[0]
    idx = np.intersect1d(idx_ra, idx_dec)
    return Q[idx], ra[idx], dec[idx]
