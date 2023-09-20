# -*- coding: utf-8 -*-
"""


@author: ia19
"""

import numpy as np
from setup_model import i,s,d,lim,r,p

def goveqs_basis2(t, in_vec, i, s, M, agg, sel, r, p):
    
    return 

in_vec = {}



out = np.zeros(len(in_vec))
invec = in_vec[:i['nstates']]

# New infections
lam = M['lam'] * invec / np.sum(invec) * (1 - p['betadec']) ** max(t - 2010, 0)

# Full model
allmat = M['lin'] + lam * M['nlin']
out[:i['nstates']] = np.dot(allmat, invec)

# Mortality
morts = M['mort'] * invec
out[:i['nstates']] -= np.sum(morts, axis=1)

allmorts = np.sum(morts)
births = p['birth'] * allmorts
out[i['U']['dom']] += births

