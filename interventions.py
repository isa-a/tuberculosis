# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 22:26:07 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import solve_ivp
from obj import get_objective
from allocate import allocate_parameters
from make_model import make_model
from gov_eqs_scaleup import goveqs_scaleup
from setup_model import i,s,r,p,agg,sel,ref,xi,prm,gps_born,likelihood
from scipy.integrate import odeint
from goveqs_basis import goveqs_basis2
import copy



# Assuming you've already written the Python equivalents of your functions
# like get_objective2, allocate_parameters, make_model, goveqs_basis2, and goveqs_scaleup

# Load data
# Assuming calibration_res.mat and Model_setup.mat have been converted to Python-friendly formats
xsto = np.load('xsto.npy')

# Function definition
obj = lambda x: get_objective(x, ref, prm, gps_born,likelihood)

# Selection of samples
ix0 = round(len(xsto) / 2)
dx = round(len(xsto) / 2 / 150)
xs = xsto[ix0::dx, :]

mk = round(len(xs) / 25)
for ii in range(len(xs)):
    if ii % mk == 0:
        print(f"{ii/mk:.5g}", end=" ")

    xx = xs[ii, :]
    out, aux = obj(xx)
    init = aux['soln'][-1, :]

    p0, r0 = allocate_parameters(xx, p, r, xi)
    M0 = make_model(p0, r0, i, s, gps_born)
    
    p1, r1 = copy.deepcopy(p0), copy.deepcopy(r0)
    p1['migrTPT'] = 1
    M1 = make_model(p0, r0, i, s, gps_born)
    
    p2, r2 = copy.deepcopy(p0), copy.deepcopy(r0)
    p2['migrTPT'] = 1
    r2['ACF'] = [0.69 * x for x in [1, 1]]
    M2 = make_model(p0, r0, i, s, gps_born)
    
    # Models list
    models = [M0, M1, M2]
    
    incsto = np.zeros((14, 167, 5))
    mrtsto = np.zeros((14, 167, 5))

    for mi, model in enumerate(models):
        def geq(t, in_):
            return goveqs_basis2(t, in_, i, s, M0, agg, sel, r, p).flatten()
        soln0 = odeint(geq, init, np.arange(2022, 2037), tfirst=True)

        sdiff = np.diff(soln0, axis=0)
        incsto[:, ii, mi] = sdiff[:,i["aux"]["inc"][0]] * 1e5
        mrtsto[:, ii, mi] = sdiff[ :,-1] * 1e5

        # ... similar code for other calculations ...

print()
