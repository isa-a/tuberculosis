# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:05:31 2023

@author: ia19
"""


born = ['dom','for']
states = ['U','Lf','Ls','I','R']
nstates = len(states)
i = []
s = []
d = []
lim = 0
fnames = []



def get_addresses(groups, i=None, s=None, d=None, lim=0):
    if i is None:
        i = {}
    if s is None:
        s = {}
    if d is None:
        d = {}

    for ig in range(len(groups)):
        gp = groups[ig]
        for ig2 in range(len(gp)):
            if gp[ig2] not in s:
                s[gp[ig2]] = []

    if len(groups) == 1:
        gp1 = groups[0]
        for ig1 in range(len(gp1)):
            lim += 1
            i[(gp1[ig1],)] = lim
            s[gp1[ig1]].append(lim)
            d[lim] = [gp1[ig1]]

    if len(groups) == 2:
        gp1 = groups[0]
        gp2 = groups[1]
        for ig1 in range(len(gp1)):
            for ig2 in range(len(gp2)):
                lim += 1
                i[(gp1[ig1], gp2[ig2])] = lim
                s[gp1[ig1]].append(lim)
                s[gp2[ig2]].append(lim)
                d[lim] = [gp1[ig1], ' ', gp2[ig2]]

    i['nstates'] = lim
    return i, s, d, lim

# Example usage:
states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R']
gps_born = ['dom', 'for']

i, s, d, lim = get_addresses([states, gps_born])

