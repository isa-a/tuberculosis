# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 00:44:17 2023

@author: ISA
"""

#SIR model implemented through linear algebra

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



# -- Natural history parameters -------------------------------------------
progression  = 0.0826
LTBI_stabil  = 0.872
reactivation = 0.0006
Tx           = 2
self_cure    = 1/6
muTB         = 1/6
prop_imm     = 0.8
# Interventions 
migrTPT      = 0
TPTeff       = 0.6                                                     

    
# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
#J0 = I0
Lf0, Ls0 = 0, 0
u = LTBI_stabil
v = progression
w = reactivation
gamma = self_cure
beta = 13
U = U0
Lf = Lf0
Ls=Ls0
I=I0
R=R0
t= np.linspace(0,100,100+1)

def get_addresses():
    # Create the matrix
    matrix = np.zeros((5, 5))
    
    # Define the row and column labels
    row_labels = ['U', 'Lf', 'Ls', 'I', 'R']
    col_labels = ['U', 'Lf', 'Ls', 'I', 'R']
    
    # Create a dictionary to map labels to indices
    index_map = {(row_labels[i], col_labels[j]): (i, j) for i in range(5) for j in range(5)}
    
    return index_map


def make_model(y, t, N, beta, gamma, u, v, w):
        ####create matrix
    U,Lf,Ls,I,R = y
    #create state vector
    state_vec = np.array([U, Lf, Ls, I, R])
    
    #create matrix full of zeros
    zero_mat = np.zeros((5,5))
    
    zero_mat[index_map[('U', 'Lf')]] = 1
    zero_mat[index_map[('Ls', 'I')]] = 2
    zero_mat[index_map[('R', 'R')]] = 3

    
    

    return

#create all matrices in here 
def make(y, t, N, beta, gamma, u, v, w):
    ####create matrix
    U,Lf,Ls,I,R = y
    #create state vector
    state_vec = np.array([U, Lf, Ls, I, R])
    
    #create matrix full of zeros
    zero_mat = np.zeros((5,5))
    
    #put the rates in correct positions
    #positions are based on a flattend matrix
    zero_mat.put([6], -u-v)
    zero_mat.put([11], u)
    zero_mat.put([12], -w)
    zero_mat.put([16], v)
    zero_mat.put([17], w)
    zero_mat.put([18], -gamma)
    zero_mat.put([23], gamma)
    
    #rename to linear component
    linearmatrix = zero_mat
    
    #create matrix of zeros for nonlinear
    #put 1s in places where lambda will be
    lam_zeros_mat = np.zeros((5,5))
    lam_zeros_mat.put([0], -1)
    lam_zeros_mat.put([5], 1)
    
    #create lambda vector, 1x5 shape
    #place beta in positions that will match up with infectious states
    lam_vector = np.zeros((5))
    lam_vector.put([3], beta)
    #get scalar lambda value by multiplying with state vector
    lamda = np.dot(lam_vector, state_vec)
    
    #add lambda to nonlinear matrix
    nonlinearmatrix = lamda * lam_zeros_mat
    
    #both nonlinear and linear components combined
    #combinedmatrices = nonlinearmatrix + linearmatrix
    
    #combined matrix multiplied by state vector to give 5x1 vector
    #solver_feed = np.dot(combinedmatrices, state_vec)
    
    #convert to tuple
    #solver_feed = solver_feed.tolist()    
    #solver_feed = tuple(solver_feed)
    
    return linearmatrix, nonlinearmatrix

    

def model_spec(y, t, N, beta, gamma, u, v, w):
    ####create matrix
    U,Lf,Ls,I,R = y
    #create state vector
    state_vec = np.array([U, Lf, Ls, I, R])
    
    #create matrix full of zeros
    zero_mat = np.zeros((5,5))    
    
    
    #get_addresses replication...will realign later
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
    
    #put the rates in correct positions
    #positions are based on a flattend matrix
    zero_mat.put([6], -u-v)
    zero_mat.put([11], u)
    zero_mat.put([12], -w)
    zero_mat.put([16], v)
    zero_mat.put([17], w)
    zero_mat.put([18], -gamma)
    zero_mat.put([23], gamma)
    
    #rename to linear component
    linearmatrix = zero_mat
    
    #create matrix of zeros for nonlinear
    #put 1s in places where lambda will be
    lam_zeros_mat = np.zeros((5,5))
    lam_zeros_mat.put([0], -1)
    lam_zeros_mat.put([5], 1)
    
    #create lambda vector, 1x5 shape
    #place beta in positions that will match up with infectious states
    lam_vector = np.zeros((5))
    lam_vector.put([3], beta)
    #get scalar lambda value by multiplying with state vector
    lamda = np.dot(lam_vector, state_vec)
    
    #add lambda to nonlinear matrix
    nonlinearmatrix = lamda * lam_zeros_mat
    
    #both nonlinear and linear components combined
    combinedmatrices = nonlinearmatrix + linearmatrix
    
    #combined matrix multiplied by state vector to give 5x1 vector
    solver_feed = np.dot(combinedmatrices, state_vec)
    
    #convert to tuple
    solver_feed = solver_feed.tolist()    
    solver_feed = tuple(solver_feed)

    return solver_feed

solve = odeint(model_spec, (U0, Lf0, Ls0, I0, R0), t, args=(N, beta, gamma, u, v, w))
U, Lf, Ls, I, R = solve.T

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, U*100000, 'black', alpha=1, lw=2, label='uninfected')
#ax.plot(t, Lf*100000, 'black', alpha=1, lw=2, label='latent fast')
#ax.plot(t, Ls*100000, 'purple', alpha=1, lw=2, label='latent slow')
ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
ax.plot(t, R*100000, 'blue', alpha=1, lw=2, label='recovered')
ax.set_xlabel('Time')
ax.set_ylabel('Number')
#ax.set_xlim(2019, 2030)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#plt.title("Incidence")
plt.show()


# Create the matrix
matrix = np.zeros((5, 5))

# Define the row and column labels
row_labels = ['U', 'Lf', 'Ls', 'I', 'R']
col_labels = ['U', 'Lf', 'Ls', 'I', 'R']

# Create a dictionary to map labels to indices
index_map = {(row_labels[i], col_labels[j]): (i, j) for i in range(5) for j in range(5)}

# Assign values to the matrix
matrix[index_map[('U', 'Lf')]] = 1
matrix[index_map[('Ls', 'I')]] = 2
matrix[index_map[('R', 'R')]] = 3

# Print the matrix
print(matrix)
