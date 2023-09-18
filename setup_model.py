# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 18:54:00 2023

@author: ISA
"""
from get_addresses import get_addresses

#combining states and born will be used as groups in get_addresses
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R'] # states
gps_born = ['dom', 'for'] # where they are born

i, s, d, lim = get_addresses([states, gps_born])

