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


gps, lhd = gps_born,likelihood

# Objectives
# Objectives
obj = lambda x: get_objective(x, ref, prm, gps, lhd)
nobj = lambda x: -obj(x)

# Number of samples
nsam = int(100)
mk = int(nsam / 25)

# Set the number of samples
nsam = 100

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
obj = lambda x: get_objective(x, ref, prm, gps_born,likelihood_function)[0]

cov0 = np.eye(len(x0))
# Perform MCMC
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, nsam, 0.1, cov0)


# Find the parameter set with the maximum log-posterior density
inds = np.where(outsto == np.max(outsto))[0]
x0 = xsto[inds[0], :]

# Run MCMC again with updated cov0 (without blockinds or fixinds)
cov0 = np.cov(xsto.T)
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, 100, 1, cov0)


xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, int(100), 1, None, None, None, True)

# Find indices of max outsto
inds = np.where(outsto == np.max(outsto))[0]
x0 = xsto[:, inds[0]]

# Evaluate objective with x0
out, aux = obj(x0)
sfin = aux['soln'][-1, :]
result = np.sum(sfin[np.intersect1d(sel['for'], [sel['Lf'], sel['Ls']])]) / np.sum(sfin[sel['for']])

# Calculate covariance
cov0 = np.cov(xsto)

# Perform MCMC again
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, int(1e4), 1, None, None, cov0, True)


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
    out, aux = obj(x)  # Assuming obj is your objective function
    # Access the values in aux using keys
    sim[ii, 0] = aux['incd'][0]  # Change this to the correct key for incd
    sim[ii, 1] = aux['mort']
    sim[ii, 2] = aux['p_migrTB']
    sim[ii, 3] = aux['p_migrpopn']
    sim[ii, 4] = aux['p_LTBI']

print('\n')