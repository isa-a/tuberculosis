# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 00:44:17 2023

@author: ISA
"""

#SIR model implemented through linear algebra

import numpy as np


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
#free params
params = ['beta', 'gamma']
#set param bounds for distribution
beta_bounds = [0, 50]
gamma_bounds = [0, 40]


####create matrix

u = LTBI_stabil
v = progression
w = reactivation
gamma = self_cure
beta = 0.4




#def matrices():
zero_mat = np.zeros((5,5))

addressU = zero_mat[0,]
addressLf = zero_mat[1,]
addressLs = zero_mat[2,]
addressI = zero_mat[3,]
addressR = zero_mat[4,]

colU = zero_mat[:,0]
colLf = zero_mat[:,1]
colLs = zero_mat[:,2]
colI = zero_mat[:,3]
colR = zero_mat[:,4]


zero_mat.put([6], -u-v)
zero_mat.put([11], u)
zero_mat.put([12], -w)
zero_mat.put([16], v)
zero_mat.put([17], w)
zero_mat.put([18], -gamma)
zero_mat.put([23], gamma)

zero_mat


lam_zeros_mat = np.zeros((5,5))
lam_zeros_mat.put([0], -1)
lam_zeros_mat.put([5], 1)


def zeroes(lamda, u, v, w, gamma):
    
    zero_mat = np.array([[-lamda, 0., 0., 0., 0.],
                         [lamda, -(u+v), 0., 0., 0.],
                         [0.,   u,      -w, 0., 0.],
                         [0.,   v,   w,     -gamma, 0.],
                         [0., 0.,   0.,      gamma, 0.]])
    
    state_vec = np.array(['U', 'Lf', 'Ls', 'I', 'R']) #needs to be changed
    
    dU = zero_mat[0,0] * state_vec[0]
    dlf = zero_mat[1,0] * state_vec[0] + zero_mat[1,1] * state_vec[1]
    dls = zero_mat[2,1] * state_vec[1] + zero_mat[2,2] * state_vec[2]
    dI = zero_mat[3,1] * state_vec[1] + zero_mat[3,2] * state_vec[2] + zero_mat[3,3] * state_vec[3]
    dR = zero_mat[4,3] * state_vec[3]
    
    #state_vec = np.array([[U, Lf, Ls, I, R]])
    return(zero_mat)

state_vec = np.array(['U', 'Lf', 'Ls', 'I', 'R'])

zeros = np.zeros((5,5))
params = ['-lambda', 'lambda', '-u-v', 'u', '-w', 'v', 'w','-gamma', 'gamma']
pos = [(0,0),(1,0),(1,1),(2,1),(2,2),(3,1),(3,2),(3,3),(4,3)]
rows, cols = zip(*pos)
zeros[rows, cols] = params


zeros = np.zeros((5, 5)).astype(str)
params = ['-lambda', 'lambda', '-u-v', 'u', '-w', 'v', 'w','-gamma', 'gamma']
posx = [0, 1, 1, 2, 2, 3, 3, 3, 4]
posy = [0, 0, 1, 1, 2, 1, 2, 3, 3]
zeros[posx, posy] = params

zeros[0,0]
#rates_matrix = np.array([[0,0,0,0,0], [0,0,0,0,0]])

print(np.dot(zeros,state_vec))




def zeroes(lamda, u, v, w, gamma):
    
    zero_mat = np.array([[-lamda, 0., 0., 0., 0.],
                         [lamda, -(u+v), 0., 0., 0.],
                         [0.,   u,      -w, 0., 0.],
                         [0.,   v,   w,     -gamma, 0.],
                         [0., 0.,   0.,      gamma, 0.]])
    
    state_vec = np.array(['U', 'Lf', 'Ls', 'I', 'R']) #needs to be changed
    
    dU = zero_mat[0,0] * state_vec[0]
    dlf = zero_mat[1,0] * state_vec[0] + zero_mat[1,1] * state_vec[1]
    dls = zero_mat[2,1] * state_vec[1] + zero_mat[2,2] * state_vec[2]
    dI = zero_mat[3,1] * state_vec[1] + zero_mat[3,2] * state_vec[2] + zero_mat[3,3] * state_vec[3]
    dR = zero_mat[4,3] * state_vec[3]
    
    #state_vec = np.array([[U, Lf, Ls, I, R]])
    return(zero_mat)

I0, R0 = 0.001, 0
lamda = beta/I0

zero_mat =     np.array([[lamda, 0., 0., 0., 0.],
                         [lamda, -(u+v), 0., 0., 0.],
                         [0., u, -w, 0., 0.],
                         [0., v, w, -gamma, 0.],
                         [0., 0., 0., gamma, 0.]])
    
state_vec = np.array(['U', 'Lf', 'Ls', 'I', 'R'])


lamdamat = np.array([0, 0, 0, beta, 0])


lamdazeros =   np.array([[-1, 0., 0., 0., 0.],
                         [1, 0, 0., 0., 0.],
                         [0., 0, 0, 0., 0.],
                         [0., 0, 0, 0, 0.],
                         [0., 0., 0., 0, 0.]])
    
dU = zero_mat[0,0] * state_vec[0]
dlf = zero_mat[1,0] * state_vec[0] + zero_mat[1,1] * state_vec[1]
dls = zero_mat[2,1] * state_vec[1] + zero_mat[2,2] * state_vec[2]
dI = zero_mat[3,1] * state_vec[1] + zero_mat[3,2] * state_vec[2] + zero_mat[3,3] * state_vec[3]
dR = zero_mat[4,3] * state_vec[3]


import numpy as np
from scipy.sparse import diags, csr_matrix

def make_model(p, r, i, s, gps):
    m = np.zeros((i['nstates'], i['nstates']))

    for ib in range(len(gps['born'])):
        born = gps['born'][ib]
        geti = lambda st: i[st][born]

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

        # Progression from 'fast' latent
        source  = Lf
        destin  = I
        rate    = r['progression']
        m[destin, source] = m[destin, source] + rate

        source  = Pf
        destin  = I2
        rate    = r['progression'] * (1 - p['TPTeff'])
        m[destin, source] = m[destin, source] + rate

        # Stabilisation of 'fast' to 'slow' latent
        source = Lf
        destin = Ls
        rate   = r['LTBI_stabil']
        m[destin, source] = m[destin, source] + rate

        source = Pf
        destin = Ps
        rate   = r['LTBI_stabil']
        m[destin, source] = m[destin, source] + rate

        # Reactivation of 'slow' latent
        source  = Ls
        destin  = I
        rate    = r['reactivation']
        m[destin, source] = m[destin, source] + rate

        source  = Ps
        destin  = I
        rate    = r['reactivation'] * (1 - p['TPTeff'])
        m[destin, source] = m[destin, source] + rate

        # Initiation of treatment
        source  = I
        destins = [Tx, Rhi]
        rates   = [r['gamma'], r['self_cure']]
        m[destins, source] = m[destins, source] + rates

        source  = I2
        destins = [Tx, Rhi]
        rates   = [r['gamma'], r['self_cure']]
        m[destins, source] = m[destins, source] + rates

        # Treatment completion or interruption
        source  = Tx
        destins = [Rlo, Rhi]
        rates   = [r['Tx'], r['default']]
        m[destins, source] = m[destins, source] + rates

        # Relapse
        sources = [Rlo, Rhi, R]
        destin  = I2
        rates   = r['relapse']
        m[destin, sources] = m[destin, sources] + rates

        # Stabilisation of relapse risk
        sources = [Rlo, Rhi]
        destin  = R
        rates   = 0.5
        m[destin, sources] = m[destin, sources]




    m = np.zeros((i['nstates'], i['nstates']))
    for ib in range(len(gps['born'])):
        born = gps['born'][ib]
        susinds = np.intersect1d([s['U'], s['Lf'], s['Ls'], s['Rlo'], s['Rhi'], s['R']], s[born])
        m[i['Lf'][born], susinds] = 1
    imminds = [s['Lf'], s['Ls'], s['Rlo'], s['Rhi'], s['R']]
    m[:, imminds] = m[:, imminds] * (1 - p['imm'])
    M['nlin'] = csr_matrix(m - np.diag(np.sum(m, axis=1)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize_scalar

# Natural history parameters
r = {'progression': 0.0826,
     'LTBI_stabil': 0.872,
     'reactivation': 0.0006,
     'self_cure': 1/6,
     'relapse': 0.003,
     'mu': 1/66,  # natural mortality
     'muTB': 1/6,  # TB related mortality
     }
p = {'imm': 0.8}  # Reduced susceptibility conferred by previous infection
prm = {'p': p, 'r': r}

x = np.array([10, 0.5])

r['beta'] = x[0]
r['gamma'] = x[1]

# --- Solve the model to equilibrium
init = np.zeros(6)
seed = 1e-6
init[0] = (1 - seed)
init[3] = seed

def geq(y, t, r, p):
    S, L, I, A, R, M = y
    dSdt = -r['beta'] * p['imm'] * S * (I + r['gamma'] * A)
    dLdt = r['beta'] * p['imm'] * S * (I + r['gamma'] * A) - r['progression'] * L - r['mu'] * L
    dIdt = r['progression'] * L - r['LTBI_stabil'] * I - r['reactivation'] * I - r['muTB'] * I
    dAdt = r['LTBI_stabil'] * I - r['self_cure'] * A - r['relapse'] * A - r['muTB'] * A
    dRdt = r['self_cure'] * A
    dMdt = r['mu'] * (S + L) + r['muTB'] * (I + A)
    return [dSdt, dLdt, dIdt, dAdt, dRdt, dMdt]

t_eval = np.linspace(0, 500, 501)
soln0 = odeint(geq, init, t_eval, args=(r, p), rtol=1e-12, atol=1e-12)

# --- Find the outputs
prev = soln0[-1, 2] * 1e5
inc = (soln0[-1, 4] - soln0[-2, 4]) * 1e5

import matplotlib.pyplot as plt
plt.plot(t_eval, soln0[:, 3] * 1e5)
plt.show()

data = {'prev': 212, 'inc': 257}
