# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born
from make_model import make_model
gps = gps_born
def goveqs_basis2(t, insert, i, s, M, agg, sel, r, p):
    
    t = 2022

    insert = np.zeros(i['nx']) # delete when finished
    seed = 1e-5
    insert[i[('U', 'dom')]] = 1 - seed
    insert[i[('I', 'dom')]] = seed


    out = np.zeros((len(insert), 1))
    invec = insert[:i['nstates']]

    # new infections
    lam = make_model(p, r, i, s, gps)['lam'] * invec / np.sum(invec) * (1 - 0.14) ** np.maximum((t - 2010), 0)


    allmat = make_model(p, r, i, s, gps)['lin'] + (lam * make_model(p, r, i, s, gps)['nlin'])
    out[:i['nstates']] = np.dot(allmat, invec.reshape(-1, 1))

    # mortality
    morts = make_model(p, r, i, s, gps)['mort'] * invec.reshape(-1, 1)
    out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)

    allmorts = np.sum(morts)
    births = 0.72 * allmorts
    out[i[('U', 'dom')]] += births

    vec = invec[s['dom']]
    vec[1:3] = vec[1:3] * 11 * (1 - p['migrTPT'])
    vec[3:5] = vec[3:5] * 11 * p['migrTPT']
    vec = vec / sum(vec) * (1 - 0.72) * allmorts
    out = out.reshape(-1)
    out[s['for']] += vec


    # aux

    return out











