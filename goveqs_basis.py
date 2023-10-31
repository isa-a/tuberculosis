# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np

def goveqs_basis2(t, insert, i, s, M, agg, sel, r, p):
        
    # initialise out vector used in odeint
    out = np.zeros((len(insert), 1))
    # select just the compartments, no aggregators or selectors
    invec = insert[:i['nstates']]
    invec = invec.reshape(-1,1)
    
    # lambda with declining beta process over time
    # will be scalar
    lam = M['lam'] @ (invec) / np.sum(invec) * (1 - p['betadec']) ** np.maximum((t - 2010), 0)
    
    # full specification of the model
    # allmat is 22x22, invec.T is 22x1
    allmat = M['lin'] + (lam * M['nlin'])
    # again just the compartments of the model
    # full model times state vec gives 22x1
    out[:i['nstates']] = np.dot(allmat, invec)
    
    # mortality
    morts = M['mort'] * invec
    # substract mortality from the states
    out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)
    
    # and births into UK pop
    dom_morts = np.sum(morts[s['dom'], :])
    out[i[('U', 'dom')]] += dom_morts
    
    
    # migration out of the UK
    out[s['migr']] = out[s['migr']] - r['migr'] * invec[s['migr']] / np.sum(invec[s['migr']])
    
    # migration in
    inmigr = np.sum(morts[s['migr'], :]) + r['migr']
    vec = np.array([
        1 - p['LTBI_in_migr'],
        (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.02,
        (1 - p['migrTPT']) * p['LTBI_in_migr'] * 0.98,
        p['migrTPT'] * p['LTBI_in_migr'] * 0.02,
        p['migrTPT'] * p['LTBI_in_migr'] * 0.98
    ])
    out[s['migrstates']] = out[s['migrstates']] + inmigr * vec
    
    
    # auxillaries    
    out[i['aux']['inc']] = agg['inc'].dot(sel['inc'] * allmat).dot(invec)
    tmp1 = agg['incsources'].dot((sel['Lf2I'] * allmat).dot(invec))
    tmp2 = agg['incsources'].dot((sel['Pf2I'] * allmat).dot(invec))
    tmp3 = agg['incsources'].dot((sel['Ls2I'] * allmat).dot(invec))
    tmp4 = agg['incsources'].dot((sel['Ps2I'] * allmat).dot(invec))
    tmp5 = agg['incsources'].dot((sel['R2I'] * allmat).dot(invec))
    out[i['aux']['incsources']] = np.concatenate([tmp1, tmp2, tmp3, tmp4, tmp5])
    out[i['aux']['mort']] = np.sum(morts[:, 1])
    out[i['aux']['nTPT']] = np.sum((sel['nTPT'] * allmat).dot(invec))
    
    return out











