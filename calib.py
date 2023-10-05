# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:58:59 2023

@author: ISA
"""

import numpy as np
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
from setup_model import i,s,lim,r,p,agg,sel,ref,xi,prm,gps_born,likelihood
from obj import get_objective
from pyDOE2 import lhs
from mcmc2 import MCMC_adaptive



# Objectives
# Objectives
obj = lambda x: get_objective(x, ref, prm, gps_born, likelihood)[0]
nobj = lambda x: -obj(x)

# Number of samples
nsam = int(1000)
mk = int(nsam / 25)

# Extract parameter names and bounds
param_names = list(prm['bounds'].keys())
param_min_values = [prm['bounds'][param][0] for param in param_names]
param_max_values = [prm['bounds'][param][1] for param in param_names]

# Generate Latin Hypercube samples
lhs_samples = lhs(len(param_names), samples=nsam)

# Scale and shift the samples to match the parameter bounds
xsam = np.zeros((nsam, len(param_names)))
for i in range(len(param_names)):
    xsam[:, i] = param_min_values[i] + lhs_samples[:, i] * (param_max_values[i] - param_min_values[i])
outs = np.zeros(nsam)

for ii in range(nsam):
    if ii % mk == 0:
        print(f'{ii / mk:.5g} ', end='')
    outs[ii] = obj(xsam[ii])
    
    
# Order by fit
mat = np.vstack((outs, np.arange(1, nsam + 1))).T
mat = mat[np.argsort(mat[:, 0])[::-1], :]
ord = mat[:, 1].astype(int)
xord = xsam[ord - 1, :]

# Initial optimization
from scipy.optimize import fmin

x0 = fmin(nobj, xord[0, :], disp=1)
x0=[  21.2257  ,  0.1499 ,   5.0908  ,  0.7880  , 21.4930]
obj = lambda x: get_objective(x, ref, prm, gps_born,likelihood)[0]
obj2 = lambda x: get_objective(x, ref, prm, gps_born,likelihood)

cov0 = np.eye(len(x0))
# Perform MCMC
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, 1000, 0.1, cov0)


# Find the parameter set with the maximum log-posterior density
inds = np.where(outsto == np.max(outsto))[0]
x0 = xsto[inds[0], :]

# Run MCMC again with updated cov0 (without blockinds or fixinds)
cov0 = np.cov(xsto.T)
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, 10000, 0.1, cov0)

#################################################################
import matplotlib.pyplot as plt

def plot_trace(xsto):
    plt.figure(figsize=(10, 6))
    plt.plot(xsto)
    plt.title("Trace plot")
    plt.xlabel("Iteration")
    plt.ylabel("Sampled value of $a$")
    plt.grid(True)
    plt.show()

# Assuming you've already run your MCMC:
# xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, 10000, 0.1, cov0)

# Plot the trace
plot_trace(xsto[:, 0])
plot_trace(xsto[:, 1])
plot_trace(xsto[:, 2])
plot_trace(xsto[:, 3])
plot_trace(xsto[:, 4])

###################################################################

# Evaluate objective with x0
out, aux = obj2(x0)
sfin = aux['soln'][-1, :]
result = np.sum(sfin[np.intersect1d(s['for'], [s['Lf'], s['Ls']])]) / np.sum(sfin[s['for']])


# Assuming you have already run the MCMC and have xsto as the parameter samples
nx = 50
ix0 = xsto.shape[0] // 2
dx = xsto.shape[0] // (2 * nx)
xs = xsto[ix0::dx, :]

mk = xs.shape[0] // 24
sim = np.zeros((xs.shape[0], 9))  # Assuming there are 9 variables in aux

for ii in range(xs.shape[0]):
    if ii % mk == 0:
        print("{:.5g} ".format(ii // mk), end="")

    x = xs[ii, :]
    out, aux = obj2(x)  # Assuming obj is your objective function
    # Access the values in aux using keys
    sim[ii, 0] = aux['incd'][0]  # Change this to the correct key for incd
    sim[ii, 1] = aux['mort']
    sim[ii, 2] = aux['p_migrTB']
    sim[ii, 3] = aux['p_migrpopn']
    sim[ii, 4] = aux['p_LTBI']

print('\n')