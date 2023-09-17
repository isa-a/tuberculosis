# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 22:20:56 2023

@author: ISA
"""

from get_addresses import i,s,d,lim
import numpy as np
from scipy.sparse import dia_matrix, csr_matrix, diags

states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R']
gps_born = ['dom', 'for']
gamma = 5


def get_states_for_born(i, born):
    state_values = {} # dict to store state values based on where they're born
    for state in states: # iterate over elements in this dict of states
        key = (state, born) # assign a key to each combo of state and born, to use as lookup in i
        state_values[state] = i.get(key) # for each state, use the key to look it up in i 
    return state_values



def make_model():
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
        rate = 0.0826  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Pf
        destin = I2
        rate = 0.0826 * (1 - 0.6)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Stabilization of 'fast' to 'slow' latent
        source = Lf
        destin = Ls
        rate = 0.872  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Pf
        destin = Ps
        rate = 0.872  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Reactivation of 'slow' latent
        source = Ls
        destin = I
        rate = 0.0006  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Ps
        destin = I
        rate = 0.0006*(1-0.6)  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Initiation of treatment
        source = I
        destins = [Tx, Rhi]
        rates = [gamma,1/6]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        source = I2
        destins = [Tx, Rhi]
        rates = [gamma,1/6]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        # Treatment completion or interruption
        source = Tx
        destins = [Rlo, Rhi]
        rates = [2, 0.01]  # Replace with random numbers
        m[destins, source] = m[destins, source] + rates
    
        # Relapse
        sources = [Rlo, Rhi, R]
        destin = I2
        rates = [0.032, 0.14, 0.0015]  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rates
    
        # Stabilization of relapse risk
        sources = [Rlo, Rhi]
        destin = R
        rates = 0.5  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rates
    
        # Initiation of TPT
        source = Lf
        destin = Pf
        rate = 0.001  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        source = Ls
        destin = Ps
        rate = 0.001  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
    
        # Case-finding
        sources = [I, I2]
        destin = Tx
        rate = 0.001  # Replace with a random number
        m[destin, sources] = m[destin, sources] + rate
    
        source = I2
        destin = Tx
        rate = 0.001  # Replace with a random number
        m[destin, source] = m[destin, source] + rate
        
    # ~~~~~~~~~~~~~~~~ LINEAR COMPONENT
    col_sums = np.sum(make_model(), axis=0)  # Calculate column sums
    modified_diagonal = make_model().diagonal() + col_sums  # Subtract column sums from diagonal
    
    # Create a sparse diagonal matrix with the modified diagonal elements
    offsets = [0]  # Diagonal offsets
    sparse_matrix = dia_matrix((modified_diagonal, offsets), shape=(i['nstates'], i['nstates']))
    sparse_matrix = sparse_matrix.toarray()
    # Now, 'sparse_matrix' contains the desired result with non-zero elements preserved
    
    M_lin = make_model() - sparse_matrix


    return M_lin
    


# ~~~~~~~~~~~~~~~~ LINEAR COMPONENT
col_sums = np.sum(make_model(), axis=0)  # Calculate column sums
modified_diagonal = make_model().diagonal() + col_sums  # Subtract column sums from diagonal

# Create a sparse diagonal matrix with the modified diagonal elements
offsets = [0]  # Diagonal offsets
sparse_matrix = dia_matrix((modified_diagonal, offsets), shape=(i['nstates'], i['nstates']))
sparse_matrix = sparse_matrix.toarray()
# Now, 'sparse_matrix' contains the desired result with non-zero elements preserved

M_lin = make_model() - sparse_matrix


# ~~~~~~~~~~~~~~~~ NON LINEAR COMPONENT

# Create an empty matrix m with the same shape as i.nstates
# Define the sizes and parameters

# Initialize the matrix m with zeros
m = np.zeros((i['nstates'], i['nstates']))

# Iterate over gps_born
for born in gps_born:
    # Find the indices of susceptible states that intersect with 'born'
    susinds = np.intersect1d([s[state] for state in ['U', 'Lf', 'Ls', 'Rlo', 'Rhi', 'R']], s[born])

    # Set the corresponding elements in m to 1
    m[i[('Lf', born)], susinds] = 1












progression  = 0.0826
LTBI_stabil  = 0.872
reactivation = 0.0006

Tx            = 2
default       = 0.01

self_cure    = 1/6
relapse      = [0.032, 0.14, 0.0015]

migrTPT      = 0                                                      
TPTeff       = 0.6
ACF          = [0, 0]
ACF2         = [0, 0]

file_path = 'matrix.txt'

# Save the matrix to the text file
np.savetxt(file_path, m, fmt='%.6f', delimiter='\t')

print(f"Matrix saved to {file_path}")


# Get the indices (coordinates) of non-zero elements in 'm'
non_zero_indices = np.transpose(np.nonzero(make_model()))

# Display the non-zero positions and their coordinates
for coord in non_zero_indices:
    print(f"Coordinate: {tuple(coord)}, Value: {make_model()[coord[0], coord[1]]}")
