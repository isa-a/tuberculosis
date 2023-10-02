# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born
from make_model import make_model
from scipy.sparse import csr_matrix


def goveqs_basis2(t, insert, i, s, M, agg, sel, r, p):
    
    #M = make_model(p, r, i, s, gps_born)
    
    # initialise out vector used in odeint
    out = np.zeros((len(insert), 1))
    # select just the compartments, no aggregators or selectors
    invec = insert[:i['nstates']]
    invec = invec.reshape(-1,1)
    
    # lambda with declining beta process over time
    # will be scalar
    lam = (M['lam'] @ invec) / np.sum(invec) * (1 - p['betadec']) ** np.maximum((t - 2010), 0)

    # full specification of the model
    # allmat is 22x22, invec.T is 22x1
    allmat = M['lin'] + (lam[0,0] * M['nlin'])
    # again just the compartments of the model
    # full model times state vec gives 22x1
    out[:i['nstates'], 0] = allmat @ invec[:,0]
    
    # mortality
    morts = M['mort'].multiply(invec)
    # substract mortality from the states
    out[:i['nstates'], 0] -= np.array(morts.sum(axis=1)).flatten()
    
    # and births
    allmorts = np.sum(morts)
    births = p['birth'] * allmorts
    # add births to uninfected compartment
    out[i[('U', 'dom')]] += births
    
    # subset of state vec selecting all domestic
    #vec = invec[s['dom']]
    vec = invec[s['dom']].flatten()
    # 1 - migrant tpt
    vec[1:3] = vec[1:3] * p['p_kLf'] * (1 - p['migrTPT'])
    # migrant tpt
    vec[3:5] = vec[3:5] * p['p_kLf'] * p['migrTPT']
    # 1- p.birth * allmorts
    vec = vec / sum(vec) * (1 - p['birth']) * allmorts
    out[s['for']] += vec.reshape(-1, 1)
    
    
    # aux
    
    out[i['aux']['inc'], 0] = (agg['inc'] @ (sel['inc'] * allmat) @ invec)[0,0]
    #out[i['aux']['inc']] = agg['inc'] @ (sel['inc'] * allmat) @ invec
    
    for j, key in enumerate(['Lf2I', 'Pf2I', 'Ls2I', 'Ps2I', 'R2I']):
        out[i['aux']['sources'][j], 0] = (csr_matrix(sel[key]).multiply(allmat) @ invec).sum()
    
    # out[i['aux']['sources'][0]] = np.sum((sel['Lf2I'] * allmat) @ invec)
    
    # out[i['aux']['sources'][1]] = np.sum((sel['Pf2I'] * allmat) @ invec)
    
    # out[i['aux']['sources'][2]] = np.sum((sel['Ls2I'] * allmat) @ invec)
    
    # out[i['aux']['sources'][3]] = np.sum((sel['Ps2I'] * allmat) @ invec)
    
    # out[i['aux']['sources'][4]] = np.sum((sel['R2I'] * allmat) @ invec)
    
    # out[-1] = np.sum(morts[:, 1])
    
    return out











