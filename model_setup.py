# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 00:44:17 2023

@author: ISA
"""

#SIR model implemented through linear algebra

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @////////////////////////////@
# @////////////////////////////@
# @/////@@@@@/////////@@@@@////@
# @///@@@@@@@@///////@@@@@@@@//@
# @@@@@@@@@@@@///////@@@@@@@@@@@
# @@@@@@@@@@@@///////@@@@@@@@@@@
# @@@@@@@@@@@@///////@@@@@@@@@@@
# @@@@@@@@@@@@///////@@@@@@@@@@@
# @@@@@@@@@@@@///////@@@@@@@@@@@
# @@@@@@@@@@@/////////@@@@@@@@@@
# @@@@@@@@///////////////@@@@@@@

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# -- Natural history parameters -------------------------------------------
progression  = 0.0826
LTBI_stabil  = 0.872
reactivation = 0.0006
Tx           = 2
self_cure    = 0.4
mu           = 1/80
muTB         = 1/6
prop_imm     = 0.8
# Interventions 
migrTPT      = 0
TPTeff       = 0.6                                                     
    
# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, U0, is susceptible to infection initially.
U0 = N - I0 - R0
#J0 = I0
Lf0, Ls0 = 0, 0
u = LTBI_stabil
v = progression
w = reactivation
gamma = self_cure
beta = 8
U = U0
Lf = Lf0
Ls=Ls0
I=I0
R=R0
t= np.linspace(0,500,500+1)


#----------------------- FUNCTIONS 1 & 2---------------------------------------
def get_addresses():
    # Create the matrix
    #n is multiplier for subpopulations
    n=1
    matrix = np.zeros((5*n,5*n))
      
    #assuming each population has 5 states
    #U, Lf, Ls, I, R
    #row_labels and col_labels represent labels for
    #rows and columns, respectively
    labels_base = ['U', 'Lf', 'Ls', 'I', 'R']
    #slice entire list
    row_labels = labels_base[:]
    col_labels = labels_base[:]
    
    #this will concatenate additional labels to the rows and cols
    #based on the number of subpops there are
    #e.g. if we have n=2, domestic and foreign born, then this will
    #create rows as U, Lf, Ls, I, R and then U_1, Lf_1, Ls_1, I_1, R_1
    for i in range((n-1)):
        row_labels += ["%s%d" % (lbl,i+1) for lbl in labels_base]
        col_labels += ["%s%d" % (lbl,i+1) for lbl in labels_base]
  
    #here I create a dictionary called index_map 
    #this iterates over all possible combinations of 
    #row and column indices using the range of whatever the matrix size is. 
    #every key in the dict is a tuple containing a pair of 
    #labels from row_labels and col_labels, while the value is a tuple of 
    #corresponding row and column indices. the dict maps the 
    #labels to their respective indices.
    index_map = {(row_labels[i], col_labels[j]): (i, j) for i in range(matrix.shape[0]) for j in range(matrix.shape[1])}
    
    return index_map


def get_lambda_addresses():
    m = 1
    vector = np.zeros((5 * m))
    vector_labels_base = ['U', 'Lf', 'Ls', 'I', 'R']
    vector_row_labels = vector_labels_base[:]
    vector_col_labels = vector_labels_base[:]
    
    for i in range((m-1)):
        vector_row_labels += ["%s%d" % (lbl,i+1) for lbl in vector_labels_base]
        vector_col_labels += ["%s%d" % (lbl,i+1) for lbl in vector_labels_base]
        
    lambda_index_map = {(vector_row_labels[i]): (i) for i in range(vector.shape[0])}
    
    return lambda_index_map


#----------------------- FUNCTION 3---------------------------------------
def make_model():
    #create matrix full of zeros
    zero_mat = np.zeros((5,5))   
    
    #use get_addresses to add rates into matrix
    #the order for the mapping is destination, source
    zero_mat[get_addresses()[('Ls', 'Lf')]] = u
    zero_mat[get_addresses()[('I', 'Lf')]] = v
    zero_mat[get_addresses()[('I', 'Ls')]] = w
    zero_mat[get_addresses()[('R', 'I')]] = gamma

    #sort out diagonal
    #function sums up columns of matrix
    def sum_diags():
        return np.sum(zero_mat, axis=0)
    
    #negative of the sum of the columns goes on the diagonal
    for index_map in zero_mat:
        np.fill_diagonal(zero_mat, -sum_diags())
    
    #call it the linear component
    linearmatrix = zero_mat
    
    #construct nonlinear component for lambda
    lam_zeros_mat = np.zeros((5,5))  
    
    #apply get addresses to add to elements where
    #the lambda values will be - add 1s 
    #call it nonlinear component
    #destination first, then source
    lam_zeros_mat[get_addresses()[('Lf', 'U')]] = 1
    nonlinearmatrix = lam_zeros_mat
    
    #sort out diagonal for nonlinear matrix
    def sum_diags_nonlin():
        return np.sum(lam_zeros_mat, axis=0)
    
    for index_map in lam_zeros_mat:
        np.fill_diagonal(lam_zeros_mat, -sum_diags_nonlin())
    
    lambda_matrix = np.zeros((5))
    lambda_matrix[get_lambda_addresses()[('I')]] = beta
    
        
    return linearmatrix, nonlinearmatrix, lambda_matrix

#----------------------- FUNCTION 4---------------------------------------


def gov_eqs(y, t, N, beta, gamma, u, v, w):
    U, Lf, Ls, I, R = y
    state_vec = np.array([U, Lf, Ls, I, R])

    mort_mat = np.zeros((5, 2))
    row_to_skip = 3

    for i in range(mort_mat.shape[0]):
        if i != row_to_skip:
            mort_mat[i, 0] = mu

    mort_mat[3, 1] = muTB

    tbdeaths = state_vec * mort_mat[:, 1]
    tbdeaths = np.sum(tbdeaths)

    births = tbdeaths / mu
    lambda_value = np.dot(make_model()[2], state_vec)

    all_mat = make_model()[0] + lambda_value * make_model()[1]

    solver_feed = np.dot(all_mat, state_vec)

    # Subtract mortality rates from the solver_feed
    solver_feed -= np.sum(mort_mat, axis=1) * state_vec
    solver_feed = np.array(solver_feed)

    # Add births to the U state
    #solver_feed[0] += births

    return solver_feed

solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), t, args=(N, beta, gamma, u, v, w))
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
ax.grid(which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#plt.title("Incidence")
plt.show()
