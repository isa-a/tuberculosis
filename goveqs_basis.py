# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born
from make_model import make_model


def goveqs_basis2(t, insert, i, s, M, agg, sel, r, p):
        
    # initialise out vector used in odeint
    out = np.zeros((len(insert), 1))
    # select just the compartments, no aggregators or selectors
    invec = insert[:i['nstates']]
    
    # lambda with declining beta process over time
    # will be scalar
    lam = make_model(p, r, i, s, gps_born)['lam'] @ (invec.reshape(-1,1)) / np.sum(invec.reshape(-1,1)) * (1 - p['betadec']) ** np.maximum((t - 2010), 0)
    
    # full specification of the model
    # allmat is 22x22, invec.T is 22x1
    allmat = make_model(p, r, i, s, gps_born)['lin'] + (lam * make_model(p, r, i, s, gps_born)['nlin'])
    # again just the compartments of the model
    # full model times state vec gives 22x1
    out[:i['nstates']] = np.dot(allmat, invec.reshape(-1, 1))
    
    # mortality
    morts = make_model(p, r, i, s, gps_born)['mort'] * invec.reshape(-1, 1)
    # substract mortality from the states
    out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)
    
    # and births
    allmorts = np.sum(morts)
    births = p['birth'] * allmorts
    # add births to uninfected compartment
    out[i[('U', 'dom')]] += births
    
    # subset of state vec selecting all domestic
    vec = invec[s['dom']]
    # 1 - migrant tpt
    vec[1:3] = vec[1:3] * p['p_kLf'] * (1 - p['migrTPT'])
    # migrant tpt
    vec[3:5] = vec[3:5] * p['p_kLf'] * p['migrTPT']
    # 1- p.birth * allmorts
    vec = vec / sum(vec) * (1 - p['birth']) * allmorts
    out[s['for']] += vec.reshape(-1, 1)
    
    
    # aux
    
    out[i['aux']['inc']] = agg['inc'] @ (sel['inc'] * allmat) @ invec.reshape(-1, 1)
    
    out[i['aux']['sources'][0]] = np.sum((sel['Lf2I'] * allmat) @ invec.reshape(-1, 1))
    
    out[i['aux']['sources'][1]] = np.sum((sel['Pf2I'] * allmat) @ invec.reshape(-1, 1))
    
    out[i['aux']['sources'][2]] = np.sum((sel['Ls2I'] * allmat) @ invec.reshape(-1, 1))
    
    out[i['aux']['sources'][3]] = np.sum((sel['Ps2I'] * allmat) @ invec.reshape(-1, 1))
    
    out[i['aux']['sources'][4]] = np.sum((sel['R2I'] * allmat) @ invec.reshape(-1, 1))
    
    out[-1] = np.sum(morts[:, 1])
    
    return out











