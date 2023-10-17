# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 16:45:39 2023

@author: ISA
"""
from setup_model import states,gps_born
import numpy as np

from scipy.sparse import csr_matrix, diags

def get_states_for_born(i, born):
    state_values = {} # dict to store state values based on where they're born
    for state in states: # iterate over elements in this dict of states
        key = (state, born) # assign a key to each combo of state and born, to use as lookup in i
        state_values[state] = i.get(key) # for each state, use the key to look it up in i 
    return state_values

def make_model2(p, r, i, s, gps):
    M = {}

    # ~~~~~~~~~~~~~~~~ LINEAR COMPONENT
    m = np.zeros((i['nstates'], i['nstates'])) # construct matrix
    for born in gps_born:
        state_values = get_states_for_born(i, born)
        
        # Access the state values for the current 'born'
        U = state_values['U']
        Lf = state_values['Lf']
        Ls = state_values['Ls']
        Pf = state_values['Pf']
        Ps = state_values['Ps']
        I = state_values['I']
        I2 = state_values['I2']
        Tx = state_values['Tx']
        Rlo = state_values['Rlo']
        Rhi = state_values['Rhi']
        R = state_values['R']
    
        # Progression from 'fast' latent
        source = Lf
        destin = I
        rate = r['progression']  # Replace with a random number
        m[destin, source] += rate
    
        source = Pf
        destin = I2
        rate = r['progression'] * (1 - p['TPTeff'])  # Replace with a random number
        m[destin, source] += rate
    
        # Stabilization of 'fast' to 'slow' latent
        source = Lf
        destin = Ls
        rate = r['LTBI_stabil']  # Replace with a random number
        m[destin, source] += rate
    
        source = Pf
        destin = Ps
        rate = r['LTBI_stabil']  # Replace with a random number
        m[destin, source] += rate
    
        # Reactivation of 'slow' latent
        source = Ls
        destin = I
        rate = r['reactivation']  # Replace with a random number
        m[destin, source] += rate
    
        source = Ps
        destin = I
        rate = r['reactivation']*(1 - p['TPTeff'])  # Replace with a random number
        m[destin, source] += rate
    
        # Initiation of treatment
        source = I
        destins = [Tx, Rhi]
        rates = [r['gamma'], r['self_cure']]  # Replace with random numbers
        m[destins, source] += rates
    
        source = I2
        destins = [Tx, Rhi]
        rates = [r['gamma'], r['self_cure']]  # Replace with random numbers
        m[destins, source] += rates
    
        # Treatment completion or interruption
        source = Tx
        destins = [Rlo, Rhi]
        rates = [r['Tx'], r['default']]  # Replace with random numbers
        m[destins, source] += rates
    
        # Relapse
        sources = [Rlo, Rhi, R]
        destin = I2
        rates = r['relapse']  # Replace with a random number
        m[destin, sources] += rates
    
        # Stabilization of relapse risk
        sources = [Rlo, Rhi]
        destin = R
        rates = 0.5  # Replace with a random number
        m[destin, sources] += rates
    
        # Initiation of TPT
        source = Lf
        destin = Pf
        rate = r['TPT'][0] # Replace with a random number
        m[destin, source] += rate
    
        source = Ls
        destin = Ps
        rate = r['TPT'][0]  # Replace with a random number
        m[destin, source] += rate
    
        # Case-finding
        sources = [I, I2]
        destin = Tx
        rate = r['ACF'][0]  # Replace with a random number
        m[destin, sources] += rate
    
        source = I2
        destin = Tx
        rate = r['ACF2'][0]   # Replace with a random number
        m[destin, source] += rate

    col_sums = np.sum(m, axis=0)  # sum up each column
    sparse_matrix = diags(col_sums, format="csr")
    M['lin'] = csr_matrix(m) - sparse_matrix

    # ~~~~~~~~~~~~~~~~ NON LINEAR COMPONENT
    m = np.zeros((i['nstates'], i['nstates']))
    for born in gps_born:
        susinds = np.intersect1d([s[state] for state in ['U', 'Ls', 'Rlo', 'Rhi', 'R']], s[born])
        m[i[('Lf', born)], susinds] = 1 

    L_and_R = [s[state] for state in ['Lf', 'Ls', 'Rlo', 'Rhi', 'R']]
    m[:, L_and_R] *= (1 - p['imm'])
    
    col_sums = np.sum(m, axis=0)
    sparse_diagonal = diags(col_sums, format="csr")
    M['nlin'] = csr_matrix(m) - sparse_diagonal

    # ~~~~~~~~~~~~~~~~ force of infection// lambda
    m = np.zeros(i['nstates'])
    m[s['everyI']] = r['beta']
    M['lam'] = diags(m, format="csr")

    # ~~~~~~~~~~~~~~~~ mortality
    m = np.zeros((i['nstates'], 2))
    m[:, 0] = 1/83
    m[s['everyI'], 1] = r['muTB']
    M['mort'] = csr_matrix(m)

    return M
