# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 22:20:56 2023

@author: ISA
"""

from get_addresses import i,s,d,lim
import numpy as np


states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R']
gps_born = ['dom', 'for']



def get_states_for_born(i, born):
    state_values = {}
    for state in states:
        key = (state, born)
        state_values[state] = i.get(key)
    return state_values


m = np.zeros((i['nstates'], i['nstates']))

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

import random


m = np.zeros((i['nstates'], i['nstates']))
    
for born in gps_born:
    state_values = get_states_for_born(i, born)
    for state in states:
        state_values = i[(state, born)]
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
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Pf
        destin = I2
        rate = random.uniform(0, 1) * (1 - 0.04)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Stabilization of 'fast' to 'slow' latent
        source = Lf
        destin = Ls
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Pf
        destin = Ps
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Reactivation of 'slow' latent
        source = Ls
        destin = I
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Ps
        destin = I
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Initiation of treatment
        source = I
        destins = [Tx, Rhi]
        rates = [random.uniform(0, 1), random.uniform(0, 1)]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        source = I2
        destins = [Tx, Rhi]
        rates = [random.uniform(0, 1), random.uniform(0, 1)]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        # Treatment completion or interruption
        source = Tx
        destins = [Rlo, Rhi]
        rates = [random.uniform(0, 1), random.uniform(0, 1)]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        # Relapse
        sources = [Rlo, Rhi, R]
        destin = I2
        rates = random.uniform(0, 1)  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rates
    
        # Stabilization of relapse risk
        sources = [Rlo, Rhi]
        destin = R
        rates = random.uniform(0, 1)  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rates
    
        # Initiation of TPT
        source = Lf
        destin = Pf
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Ls
        destin = Ps
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Case-finding
        sources = [I, I2]
        destin = Tx
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rate
    
        source = I2
        destin = Tx
        rate = random.uniform(0, 1)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    




# Assuming you have your 'm' matrix defined
# m = np.zeros((i['nstates'], i['nstates']))

# Specify the file path where you want to save the matrix
file_path = 'matrix.txt'

# Save the matrix to the text file
np.savetxt(file_path, m, fmt='%.6f', delimiter='\t')

print(f"Matrix saved to {file_path}")
