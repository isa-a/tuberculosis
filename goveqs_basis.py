# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born
from make_model import make_model
gps = gps_born

def goveqs_basis2():

    return 


    
t = 2022

insert = np.zeros(i['nx']) # delete when finished
seed = 1e-5
insert[i[('U', 'dom')]] = 1 - seed
insert[i[('I', 'dom')]] = seed


out = np.zeros((len(insert), 1))
invec = insert[:i['nstates']]

# new infections
lam = make_model(p, r, i, s, gps)['lam'] @ (invec.reshape(-1,1)) / np.sum(invec.reshape(-1,1)) * (1 - 0.1473)


allmat = make_model(p, r, i, s, gps)['lin'] + (lam * make_model(p, r, i, s, gps)['nlin'])
out[:i['nstates']] = np.dot(allmat, invec.reshape(-1, 1))

# mortality
morts = make_model(p, r, i, s, gps)['mort'] * invec.reshape(-1, 1)
out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)

# and births
allmorts = np.sum(morts)
births = 0.7268 * allmorts
out[i[('U', 'dom')]] += births

vec = invec[s['dom']]
vec[1:3] = vec[1:3] * 12.3779 * (1 - 0)
vec[3:5] = vec[3:5] * 11 * 0
vec = vec / sum(vec) * (1 - 0.72) * allmorts
out[s['for']] += vec.reshape(-1, 1)


# aux

out[i['aux']['inc']] = agg['inc'] @ (sel['inc'] * allmat) @ invec.reshape(-1, 1)

out[i['aux']['sources'][0]] = np.sum((sel['Lf2I'] * allmat) @ invec.reshape(-1, 1))

out[i['aux']['sources'][1]] = np.sum((sel['Pf2I'] * allmat) @ invec.reshape(-1, 1))

out[i['aux']['sources'][2]] = np.sum((sel['Ls2I'] * allmat) @ invec.reshape(-1, 1))

out[i['aux']['sources'][3]] = np.sum((sel['Ps2I'] * allmat) @ invec.reshape(-1, 1))

out[i['aux']['sources'][4]] = np.sum((sel['R2I'] * allmat) @ invec.reshape(-1, 1))

out[i['aux']['mort']] = np.sum(morts[:, 1])










