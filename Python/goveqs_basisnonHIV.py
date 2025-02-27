

import numpy as np
from scipy.sparse import csr_matrix
import pdb
from setup_model import *

def goveqs_basisnonHIV(t, in_vec, i, s, M, agg, sel, r, p, prm):


    out = np.zeros(i['nx'])
    invec = in_vec[: i['nstates']]
    invec = np.array(invec, dtype=float)

    # Prepare population denominators
    tmp = M['denvec'] @ invec
    den = np.sum(M['denvec'].toarray() * tmp[:, None], axis=0)
    den[den == 0] = np.inf

    try:
        lam_val = (M['lam'] @ (invec / den)) * ((1 - p['betadec']) ** max(t - 2010, 0))
        allmat = M['lin'] + lam_val[0] * M['nlin']['ch'] + lam_val[1] * M['nlin']['ad']
        out = allmat.dot(invec)
    except Exception as e:
        print("Error in computing non-HIV equations:", e)
        raise

        
    # Mortality and births
    morts = M['mort'].toarray() * invec.reshape(-1, 1)
    out[: i['nstates']] = out[: i['nstates']] - np.sum(morts, axis=1)

    dom_morts = np.sum(morts[np.array(s['dom']) - 1, :])
    out[i['U']['ch']['dom']['neg'] - 1] = out[i['U']['ch']['dom']['neg'] - 1] + dom_morts
    
        # Convert arrays to numeric types
    sel_inc = np.array(sel['inc'], dtype=float)
    allmat_f = allmat.toarray()
    invec_f = np.array(invec, dtype=float)


    # Auxiliaries
    #out[np.array(i['aux']['inc']) - 1] = agg['inc'] @ (np.multiply(sel_inc, allmat_f) @ invec_f)
# Convert the auxiliary indices to integers and shift to 0-indexing
    aux_indices = np.array(i['aux']['inc'], dtype=int) - 1
    # Force any index that's too high to equal i['nx'] - 1
    aux_indices[aux_indices >= i['nx']] = i['nx'] - 1
    # Now assign
    out[aux_indices] = agg['inc'] @ (np.multiply(sel_inc, allmat_f) @ invec_f)


    out[np.array(i['aux']['mort']) - 1] = np.sum(morts[:, 1])
    out[np.array(i['aux']['nTPT']) - 1] = (sel['nTPT'].multiply(allmat) @ invec).sum()
    out[np.array(i['aux']['ch_notifs']) - 1] = (sel['ch_notifs'].multiply(allmat) @ invec).sum()


    return out


