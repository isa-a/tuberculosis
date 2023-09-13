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

def get_the_addresses(groups, i, s, d, lim):
    
    for item in range(1, len(groups)):
        gp = groups(item)
        for item2 in range(1, len(gp)):
            if 
      
    
    
    if len(groups) == 2:
        group1 = groups[0]
        
        for item in range(1, len(group1)):            
            lim += 1
            i[group1[item - 1]] = lim  # Note: Python uses 0-based indexing
            s[group1[item - 1]] = s.get(group1[item - 1], []) + [lim]  # Initialize with an empty list if necessary
            d[lim] = [group1[item - 1]]

    return i,s,d,lim

get_the_addresses([states, born], i, s, d, lim) == i,s,d,lim



def get_addresses(groups):
    i = {}
    s = {}
    d = []
    lim = 0

    # Initiate any sets not covered by s
    if s:
        fnames = s.keys()
    else:
        fnames = set()

    for gp in groups:
        for item in gp:
            if item not in fnames:
                s[item] = []

    for group_length in range(1, len(groups) + 1):
        for combo in itertools.product(*groups[:group_length]):
            lim += 1
            if group_length == 1:
                i[combo[0]] = lim
                s[combo[0]].append(lim)
                d.append([combo[0]])
            else:
                key = '.'.join(combo)
                if key not in i:
                    i[key] = lim
                s[combo[0]].append(lim)
                s[combo[1]].append(lim)
                d.append([' '.join(combo)])

    i['nstates'] = lim
    return i, s, d

import itertools

states = ['U', 'Lf', 'Ls', 'I', 'R']
gps_born = ['dom', 'for']

i, s, d = get_addresses([states, gps_born])
d = [' '.join(item) for item in d]


print(i)
print(s)
print(d)
print(lim)




def get_addresses(groups, i, s, d, lim):
    # Rest of the code remains the same
    if not s:
        fnames = []
    else:
        fnames = s.keys()

    for gp in groups:
        for item in gp:
            if item not in fnames:
                s[item] = []

    for group_length in range(1, len(groups) + 1):
        for combo in itertools.product(*groups[:group_length]):
            lim += 1
            if group_length == 1:
                i[combo[0]] = lim
                s[combo[0]].append(lim)
                d.append([combo[0]])
            else:
                key = '.'.join(combo)
                if key not in i:
                    i[key] = lim
                s[combo[0]].append(lim)
                s[combo[1]].append(lim)
                d.append([' '.join(combo)])

    i['nstates'] = lim
    return i, s, d, lim

import itertools

states = ['U', 'Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Rlo', 'Rhi', 'R']
gps_born = ['dom', 'for']

i, s, d, lim = get_addresses([states, gps_born], {}, {}, [], 0)
d = [' '.join(item) for item in d]

print(i)
print(s)
print(d)












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



# Print the results
print("i:", result_i)
print("s:", result_s)
print("d:", result_d)
print("lim:", result_lim)

