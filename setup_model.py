# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 18:54:00 2023

@author: ISA
"""
from get_addresses import get_addresses
import numpy as np
from scipy.sparse import csr_matrix, dia_matrix



# ~~~~~~~~~~~ Natural history params

# all rates
r = {
    'progression': 0.0826,
    'LTBI_stabil': 0.872,
    'reactivation': 0.0006,
    'Tx': 2,
    'default': 0.01,
    'self_cure': 1/6,
    'relapse': [0.032, 0.14, 0.0015],
    'muTB': 1/6
    }

# all proportions
p = {
    'imm': 0.8,
    'migrTPT': 0,
    'TPTeff': 0.6
    }

sel = {}
agg = {}


r['TPT']          = [0,0]
r['ACF']          = [0, 0]
r['ACF2']         = [0, 0]


# combining states and born will be used as groups in get_addresses
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R'] # states
gps_born = ['dom', 'for'] # where they are born

i, s, d, lim = get_addresses([states, gps_born])

s['everyI'] = np.concatenate((s['I'], s['I2']))



# ~~~~~~~~~~~~ auxillary time


auxillaries = ['inc', 'sources', 'mort']
lengths = [3, 5, 1]
lim=i['nstates']
i['aux'] = {} #initialise aux in i

for ii in range(len(auxillaries)): # loop over however many auxillaries there are
    # in each iteration of the loop create a list of integers starting 
    # from lim + 1 and ending at lim + lengths[ii] + 1. Here, lim is a variable that 
    # keeps track of the last index used, and lengths[ii] is the length associated 
    # with the current auxiliary.
    inds = list(range(lim + 1, lim + lengths[ii] + 1))
    i['aux'][auxillaries[ii]] = inds
    lim = inds[-1]
i['nx'] = lim

# --- Incidence
tmp = np.zeros((3, i['nstates']))
tmp[0, s['everyI']] = 1
tmp[1, np.intersect1d(s['everyI'], s['dom'])] = 1
tmp[2, np.intersect1d(s['everyI'], s['for'])] = 1
agg['inc'] = csr_matrix(tmp).toarray()

tmp = np.zeros((i['nstates'], i['nstates']))
tmp[s['everyI'], :] = 1
sel['inc'] = tmp - np.diag(np.diag(tmp)) # so that diagonal self to self terms arent counted

# Create empty dictionaries to store selectors and aggregators

# --- Incidence

# From recent infection
tmp = np.zeros((i['nstates'], i['nstates']))
tmp[0, s['everyI']] = 1
sel['Lf2I'] = tmp - np.diag(np.diag(tmp))

# From recent infection, TPT
tmp = np.zeros((i['nstates'], i['nstates']))
tmp[1, s['everyI']] = 1
sel['Pf2I'] = tmp - np.diag(np.diag(tmp))

# From remote infection
tmp = np.zeros((i['nstates'], i['nstates']))
tmp[0, s['everyI']] = 1
sel['Ls2I'] = tmp - np.diag(np.diag(tmp))

# From remote infection, TPT
tmp = np.zeros((i['nstates'], i['nstates']))
tmp[1, s['everyI']] = 1
sel['Ps2I'] = tmp - np.diag(np.diag(tmp))

# From relapse
tmp = np.zeros((i['nstates'], i['nstates']))
tmp[0, s['everyI']] = 1
tmp[1, [s['R'][0], s['Rlo'][0], s['Rhi'][0]]] = 1
sel['R2I'] = tmp - np.diag(np.diag(tmp))




#~~~~~~~~~~~~~~~~~ parameters

free_params = ['beta', 'betadec', 'gamma', 'p_birth', 'p_kLf']
param_lengths = [1,         1,         1,       1,          1]


limit = 0
xi = {}

for ii in range(len(free_params)):
    indices = list(range(limit + 1, limit + param_lengths[ii] + 1))
    xi[free_params[ii]] = indices
    limit = indices[-1]

prm = {}
prm['bounds'] = {
    'beta': [0, 40],
    'betadec': [0, 0.15],
    'gamma': [0, 50],
    'p_birth': [0, 1],
    'p_kLf': [1, 200]
}

ref = {}


ref['i'] = i
ref['s'] = s
ref['xi'] = xi


# ~~~~~~~~~~~~~~~~~~~~~data 

data = {
    'incd2010': [14.1, 14.6, 15.1],
    'incd2020': [6.5, 7, 7.5],
    'mort': [0.28, 0.3, 0.32],
    'p_migrTB': [0.708, 0.728, 0.748],
    'p_migrpopn': [0.138, 0.168, 0.198],
    'p_LTBI': [0.15, 0.2, 0.25]
}










