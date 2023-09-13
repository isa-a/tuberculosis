# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 22:20:56 2023

@author: ISA
"""

from get_addresses import i,s,d,lim
import numpy as np

def make_model():
    
    m = np.zeros((i['nstates'],i['nstates']))
    
    for ip in range(len(gps_born)):
        born = gps_born[ip]




def make_model(p, r, i, s, gps):
    # Create an empty matrix M with the same shape as m
    M = np.zeros((i['nstates'], i['nstates']))

    for born in gps['born']:
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

