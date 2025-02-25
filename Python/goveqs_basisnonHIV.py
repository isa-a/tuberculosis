
from make_model import M
import numpy as np
from scipy.sparse import csr_matrix
import pdb
from setup_model import *

def goveqs_basisnonHIV(t, in_vec, i, s, M, agg, sel, r, p, prm):


    out = np.zeros(len(in_vec))
    invec = in_vec[: i['nstates']]

    # Prepare population denominators
    tmp = M['denvec'] @ invec
    den = np.sum(M['denvec'].toarray() * tmp[:, None], axis=0)
    den[den == 0] = np.inf

    try:
        lam_val = (M['lam'] @ (invec / den)) * ((1 - p['betadec']) ** max(t - 2010, 0))
        allmat = (M['lin'] +
                  lam_val[0] * M['nlin']['ch'] +
                  lam_val[1] * M['nlin']['ad'])
        out = allmat @ invec
    except Exception as e:
        pdb.set_trace()

    return out


# Mortality and births
morts = M['mort'].toarray() * invec.reshape(-1, 1)
out[: i['nstates']] = out[: i['nstates']] - np.sum(morts, axis=1)

dom_morts = np.sum(morts[np.array(s['dom']) - 1, :])
out[i['U']['ch']['dom']['neg'] - 1] = out[i['U']['ch']['dom']['neg'] - 1] + dom_morts

# Auxiliaries
out[np.array(i['aux']['inc']) - 1] = agg['inc'] @ (np.multiply(sel['inc'], allmat) @ invec)
out[np.array(i['aux']['mort']) - 1] = np.sum(morts[:, 1])
out[np.array(i['aux']['nTPT']) - 1] = np.sum((np.multiply(sel['nTPT'], allmat) @ invec))
out[np.array(i['aux']['ch_notifs']) - 1] = np.sum((np.multiply(sel['ch_notifs'], allmat) @ invec))
