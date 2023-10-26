# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 16:51:37 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
from make_model import make_model
from allocate import allocate_parameters
from setup_model import likelihood_function
from goveqs_basis import goveqs_basis2

def get_objective(x, ref, prm, gps, calfn):
    
    # extract dictionaries
    i = ref['i']
    s = ref['s']
    xi = ref['xi']
    p = prm['p']
    r = prm['r']
    sel = prm['sel']
    agg = prm['agg']
    
    
    # defined final parameters to test code - will be deleted after
    #x = [23.9410, 0.1473, 5.2634, 0.7268, 12.3779]
    # assign parameters to p and r dicts
    p, r = allocate_parameters(x, p, r, xi)
    
    first_values = [prm['bounds'][key][0] for key in prm['bounds']]
    second_values = [prm['bounds'][key][1] for key in prm['bounds']]
    
    tmp1 = np.array([first_values, second_values])
    tmp1 = np.vstack((tmp1, x))
    # Calculate differences between consecutive rows for columns 0, 2, and 1 (assuming 0-based indexing)
    tmp2 = np.diff(tmp1[[0, 2, 1], :], axis=0)
    
    # Check if the minimum value in tmp2 is less than 0
    cond1 = np.min(tmp2) < 0
    
    if cond1:
        out = -np.inf
        aux = np.nan
    else:
        # creates model with params
        M = make_model(p, r, i, s, gps)
        
        # initialise the 'insert' for gov eqs basis
        # (initial state vector with aux and selectors on end)
        init = np.zeros(i['nx'])
        seed = 1e-5
        init[i[('U', 'dom')]] = 1 - seed
        init[i[('I', 'dom')]] = seed
        
        def r_TPT_linear_increase(t, r_migrTPT2019):
            if t < 2015:
                return 0
            elif 2015 <= t <= 2019:
                return (r_migrTPT2019 / (2019 - 2014)) * (t - 2014)
            else:
                return r_migrTPT2019
        
        def geq(t, in_):
            r['TPT'][1] = r_TPT_linear_increase(t, 0.3)
            return goveqs_basis2(t, in_, i, s, M, agg, sel, r, p).flatten()

        
        # wrapper for gov eqs basis taking only insert and time
        # def geq(t, in_):
        #     return goveqs_basis2(t, in_, i, s, M, agg, sel, r, p).flatten()
        
        # time range for solving equation until
        t0 = np.arange(2020)
        soln0 = odeint(geq, init, t0, tfirst=True)
        #soln0 = solve_ivp(geq, (t0[0], t0[-1]), init, t_eval=t0, vectorized=True)
        #soln0=soln0.y.T

        
        # calculates the differences between rows of soln0
        dsol = np.diff(soln0, axis=0)
        # filters rows of dsol, the t0 index select rows where the corresponding 
        # time point in t0 is equal to 2010 
        # aux index is position of incidence
        incd2010 = dsol[t0[:-1] == 2010, i['aux']['inc'][0]] * 1e5
        # final position in row of dsol gives 2020, and aux index for incidence
        incd2020 = dsol[-1, i['aux']['inc'][0]] * 1e5
        # incidence in all positions of the aux for 2020
        incd = dsol[-1, i['aux']['inc']] * 1e5
        # final mortality (2020) and its index at end of aux
        mort = dsol[-1, -1] * 1e5
        p_migrTB = incd[2] / incd[0]
        
        # extract final row of soln0, which has the solutions at the last time point
        # (selecting last row 2020, and all 31 columns)
        sfin = soln0[-1, :]
        # within this subset, look for latent foreign states, to get proportion of latent
        p_LTBI = np.sum(sfin[np.intersect1d(s['for'], [s['Lf'], s['Ls']])]) / np.sum(sfin[s['for']])
        # proportion of migrants is all foreign compartments over total pop
        p_migrpopn = np.sum(sfin[s['for']]) / np.sum(sfin[:i['nstates']])
        
        # number of tpt
        n_tpt2019 = sfin[s['Pf']] + sfin[s['Ps']]
        
        if np.any(incd > 0.1):
            #out = calfn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI)
            out = likelihood_function(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, n_tpt2019)
            # add to the auxillaries with calculated values
            aux = {
                'soln': soln0,
                'incd': dsol[np.where(t0 == 2010)[0][0]:, i['aux']['inc'][0]] * 1e5,
                'incd2010': incd2010,
                'incd2020': incd2020,
                'mort': mort,
                'p_migrTB': p_migrTB,
                'p_migrpopn': p_migrpopn,
                'p_LTBI': p_LTBI,
                'n_tpt2019': n_tpt2019
            }
        else:
            out = -np.inf
            aux = np.nan
        
    return out, aux


