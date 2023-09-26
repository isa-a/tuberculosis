# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 16:51:37 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
from make_model import make_model
from allocate import allocate_parameters
from setup_model import ref, prm, sel, agg, gps_born, likelihood_function
from goveqs_basis import goveqs_basis2
from get_addresses import get_addresses

gps = gps_born
# def get_objective2(x, ref, prm, gps, calfn):

#     return out, aux
def get_objective(x, ref, prm, gps, calfn):
    i = ref['i']
    s = ref['s']
    xi = ref['xi']
    p = prm['p']
    r = prm['r']
    sel = prm['sel']
    agg = prm['agg']
    
    x = [23.9410, 0.1473, 5.2634, 0.7268, 12.3779]
    p, r = allocate_parameters(x, p, r, xi)
    
    M = make_model(p, r, i, s, gps)
    
    
    init = np.zeros(i['nx'])
    seed = 1e-5
    init[i[('U', 'dom')]] = 1 - seed
    init[i[('I', 'dom')]] = seed
    
    def geq(t, in_):
        return goveqs_basis2(t, in_, i, s, M, agg, sel, r, p).flatten()
    
    t0 = np.arange(2021)  # Adjust the time points accordingly
    soln0 = odeint(geq, init, t0, tfirst=True)
    
    
    dsol = np.diff(soln0, axis=0)
    incd2010 = dsol[t0[:-1] == 2010, i['aux']['inc'][0]] * 1e5
    incd2020 = dsol[-1, i['aux']['inc'][0]] * 1e5
    incd = dsol[-1, i['aux']['inc']] * 1e5
    mort = dsol[-1, -1] * 1e5
    p_migrTB = incd[2] / incd[0]
    
    sfin = soln0[-1, :]
    p_LTBI = np.sum(sfin[np.intersect1d(s['for'], [s['Lf'], s['Ls']])]) / np.sum(sfin[s['for']])
    
    sfin = soln0[-1, :]
    p_migrpopn = np.sum(sfin[s['for']]) / np.sum(sfin[:i['nstates']])
    
    if np.any(incd > 0.1):
        #out = calfn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI)
        out = likelihood_function(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI)
        aux = {
            'soln': soln0,
            'incd': dsol[np.where(t0 == 2010)[0][0]:, i['aux']['inc'][0]] * 1e5,
            'incd2010': incd2010,
            'incd2020': incd2020,
            'mort': mort,
            'p_migrTB': p_migrTB,
            'p_migrpopn': p_migrpopn,
            'p_LTBI': p_LTBI,
        }
    else:
        out = -np.inf
        aux = np.nan
        
        return out, aux