# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 18:54:00 2023

@author: ISA
"""
from get_addresses import get_addresses

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

r['TPT']          = [0,0]
r['ACF']          = [0, 0]
r['ACF2']         = [0, 0]


# combining states and born will be used as groups in get_addresses
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R'] # states
gps_born = ['dom', 'for'] # where they are born

i, s, d, lim = get_addresses([states, gps_born])

s['everyI'] = [s['I'], s['I2']] # all infectious


# ~~~~~~~~~~~~ auxillary time


auxillaries = ['inc', 'sources', 'mort']
lengths = [3, 5, 1]
lim=0
i['aux'] = {} #initialise aux in i

for ii in range(len(auxillaries)): # loop over however many auxillaries there are
    inds = list(range(lim + 1, lim + lengths[ii] + 1))
    i['aux'][auxillaries[ii]] = inds
    lim = inds[-1]
i['nx'] = lim





























