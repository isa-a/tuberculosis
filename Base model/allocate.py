# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 13:29:40 2023

@author: ISA
"""

from setup_model import p,r,xi



def allocate_parameters(x, p, r, xi):
    r['beta'] = x[xi['beta'][0]]
    p['betadec'] = x[xi['betadec'][0]]
    r['gamma'] = x[xi['gamma'][0]]
    p['birth'] = x[xi['p_birth'][0]]
    p['p_kLf'] = x[xi['p_kLf'][0]]
    
    return p, r
