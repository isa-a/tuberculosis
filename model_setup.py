# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 00:44:17 2023

@author: ISA
"""

#SIR model implemented through linear algebra

import numpy as np

compartments = ['U', 'Lf', 'Ls', 'I', 'R']
population = ['dom', 'for']

# -- Natural history parameters -------------------------------------------
progression  = 0.0826
LTBI_stabil  = 0.872
reactivation = 0.0006

Tx            = 2
default       = 0.01

self_cure    = 1/6
relapse      = [0.032, 0.14, 0.0015]
muTB         = 1/6
prop_imm     = 0.8

# Interventions 
migrTPT      = 0
TPTeff       = 0.6                                                     
TPT          = [0, 0]
ACF          = [0, 0]
ACF2         = [0, 0]