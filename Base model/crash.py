# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:09:30 2023

@author: ISA
"""

import numpy as np
import numpy.random as npr
from scipy.stats import multivariate_normal
from obj import get_objective
from setup_model import ref, prm, gps_born, likelihood_function


x0 = [21.2257 ,   0.1499 ,   5.0908   , 0.7880  , 21.4930]

def MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, cov0, displ):
    d = len(x0)
    b = 0.05
    sd = sigma * 2.4**2 / d

    if fixinds:
        inds = fixinds[0]
        vals = fixinds[1]
    else:
        inds = []
        vals = []

    # Checks on the initial covariance matrix
    if cov0 is None:
        cov0 = np.eye(d)
        cov0[inds, :] = 0
        cov0[:, inds] = 0
    cov0 = np.array(cov0)
    cov0[0:blockind, blockind + 1:] = 0
    cov0[blockind + 1:, 0:blockind] = 0

    # Initialize output matrices
    xsto = np.zeros((d, n))
    outsto = np.zeros(n)
    history = np.zeros((n, d + 1))  # Rows: 1:d Proposed values 4. Accept or reject

    xsto[:, 0] = x0[:]
    xbar = xsto.copy()
    FX = F(x0)
    outsto[0] = FX
    acc = 0

    if displ:
        # Add code here for any specific display options
        pass

    # --- Start the MCMC loop -------------------------------------------------
    for t in range(1, n):
        try:
            X = xsto[:, t - 1]

            # --- Make a proposal from the distribution
            Y0 = npr.multivariate_normal(X, 0.1**2 * cov0 * sigma / d)
            if t < 2 * d:
                Y = np.maximum(Y0, 0)
                Y[inds] = vals
            else:
                ind0 = 0
                ind1 = t - 1
                covmat = np.cov(xsto[:, ind0:ind1].T)
                covmat[inds, :] = 0
                covmat[:, inds] = 0
                covmat[0:blockind, blockind + 1:] = 0
                covmat[blockind + 1:, 0:blockind] = 0
                covmat = (covmat + covmat.T) / 2

                # Precaution to make sure values don't get negative
                Y = np.maximum((1 - b) * npr.multivariate_normal(X, sd * covmat) + b * Y0, 0)
                Y[inds] = vals
            history[t, 0:d] = Y

            # --- Decide whether to accept or not
            FY = F(Y)
            if (npr.rand() < np.exp(FY - FX)) and (abs(FY) < np.inf):
                # Accept
                xsel = Y[:]
                FX = FY
                acc = acc + 1
                history[t, -1] = 1
            else:
                # Reject
                xsel = xsto[:, t - 1]
            xsto[:, t] = xsel
            outsto[t] = FX
            xbar[:, t] = (xbar[:, t - 1] * (t - 1) + xsel) / t

            # Display options
            if displ and (t % round(n / 25) == 0):
                print(f'{t / n * 25:.5g} ', end='')
            if displ and (t % 200 == 0):
                # Add code here for specific display options (e.g., plot)
                pass
        except:
            # Add error handling code here
            pass

    accept_rate = acc / n
    return xsto, outsto, history, accept_rate

def obj(x):
    return get_objective(x, ref, prm, gps_born, likelihood_function)


# Set initial parameters and other inputs
x0 = np.array([0.0, 0.0, 0.0])
n_iterations = 1
initial_covariance_matrix = None
scaling_factor = 1.0
fixed_indices = None
block_indices = 0
display_options = True

xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, n_iterations, scaling_factor, fixed_indices, block_indices, initial_covariance_matrix, display_options)
