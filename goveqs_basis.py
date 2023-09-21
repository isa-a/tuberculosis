# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p
from make_model import make_model

def goveqs_basis2(t, insert, i, s, M, agg, sel, r, p):
    return 

t = 2022

insert = np.zeros(i['nx']) # delete when finished
seed = 1e-5
insert[i[('U', 'dom')]] = 1 - seed
insert[i[('I', 'dom')]] = seed


out = np.zeros((len(insert), 1))
invec = insert[:i['nstates']]

# new infections
lam = make_model()['lam'] * invec / np.sum(invec) * (1 - 0.14) ** np.maximum((t - 2010), 0)


allmat = make_model()['lin'] + (lam * make_model()['nlin'])
out[:i['nstates']] = np.dot(allmat, invec.reshape(-1, 1))


morts = make_model()['mort'] * invec.reshape(-1, 1)

out[:i['nstates']] -= np.sum(morts, axis=1).reshape(-1, 1)

allmorts = np.sum(morts)
births = 0.72 * allmorts
out[i[('U', 'dom')]] += births









