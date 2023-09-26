# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:33:40 2023

@author: ISA
"""

import numpy as np
from goveqs_basis import goveqs_basis2
from allocate import allocate_parameters
from setup_model import r,p,gps_born,i,s
from make_model import make_model



times = [2022, 2025]
t=np.arange(2022, 2037)

p['migrTPT'] = 0
p1 = p
r1 = r
p1['migrTPT'] = 1
M0 = make_model(p, r, i, s, gps_born)
M1 = make_model(p1, r1, i, s, gps_born)


def goveqs_scaleup(t, insert, i, s, M0, M1, p0, p1, times, agg, sel, r):
    #scale = min(max((t - times[0]) / (times[1] - times[0]), 0), 1)
    scale = np.minimum(1.0, np.maximum(0.0, (t - times[0]) / (times[1] - times[0])))

    
    Ms = M1.copy()
    Ms['lin'] = M0['lin'] + scale * (M1['lin'] - M0['lin'])
    
    ps = p1.copy()
    ps['migrTPT'] = p0['migrTPT'] + scale * (p1['migrTPT'] - p0['migrTPT'])
    
    out = goveqs_basis2(t, insert, i, s, Ms, agg, sel, r, ps)
    
    return out
