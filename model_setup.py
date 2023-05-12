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
# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
J0 = I0
Lf0, Ls0 = 0, 0

U = U0
Lf = Lf0
Ls=Ls0
I=I0
R=R0

#create state vector
state_vec = np.array([U, Lf, Ls, I, R])

#def matrices():
zero_mat = np.zeros((5,5))

# addressU = zero_mat[0,]
# addressLf = zero_mat[1,]
# addressLs = zero_mat[2,]
# addressI = zero_mat[3,]
# addressR = zero_mat[4,]

# colU = zero_mat[:,0]
# colLf = zero_mat[:,1]
# colLs = zero_mat[:,2]
# colI = zero_mat[:,3]
# colR = zero_mat[:,4]


zero_mat.put([6], -u-v)
zero_mat.put([11], u)
zero_mat.put([12], -w)
zero_mat.put([16], v)
zero_mat.put([17], w)
zero_mat.put([18], -gamma)
zero_mat.put([23], gamma)

linearmatrix = zero_mat


lam_zeros_mat = np.zeros((5,5))
lam_zeros_mat.put([0], -1)
lam_zeros_mat.put([5], 1)

lam_vector = np.zeros((5))
lam_vector.put([3], beta)
lamda = np.dot(lam_vector, state_vec)

nonlinearmatrix = lamda * lam_zeros_mat

combinedmatrices = nonlinearmatrix + linearmatrix

np.dot(combinedmatrices, state_vec)





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

