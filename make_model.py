# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 22:20:56 2023

@author: ISA
"""

from get_addresses import i,s,d,lim
import numpy as np

state = 'U'
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R']
gps_born = ['dom', 'for']

def make_model():
    
    m = np.zeros((i['nstates'],i['nstates']))
    
    for ip in range(len(gps_born)):
        for born in gps_born[]:
            def geti(i, st, born):
                key = (st, born)
                return i.get(key)
            
            
            U = geti(i, 'U', born)





def make_model(p, r, i, s, gps):
    # Create an empty matrix M with the same shape as m
    M = np.zeros((i['nstates'], i['nstates']))

    
        # Define a helper function to retrieve indices
        def geti(state):
            return i[state][born]

        U   = geti('U')
        Lf  = geti('Lf')
        Ls  = geti('Ls')
        Pf  = geti('Pf')
        Ps  = geti('Ps')
        I   = geti('I')
        I2  = geti('I2')
        Tx  = geti('Tx')
        Rlo = geti('Rlo')
        Rhi = geti('Rhi')
        R   = geti('R')

        # Perform operations using the retrieved indices if needed

    return M

if result is not None:
    print(result)  # This will print 1
else:
    print("Key not found")


def geti(i, st, born):
    key = (st, born)
    return i.get(key)

def make_model():
    m = np.zeros((i['nstates'], i['nstates']))
    
    for born in gps_born:
        for state in states:
            U = geti(i, state, born)
            # Do something with U
            print(U)  # Example usage

make_model()
