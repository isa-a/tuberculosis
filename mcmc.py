# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 21:21:26 2023

@author: ISA
"""

import numpy as np
import numpy.random as npr
from scipy.stats import multivariate_normal
from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born,likelihood
from obj import get_objective


# Example usage with data as lists or arrays
incd2010 = [14.1, 14.6, 15.1]
incd2020 = [6.5, 7, 7.5]
mort = [0.28, 0.3, 0.32]
p_migrTB = [0.708, 0.728, 0.748]
p_migrpopn = [0.138, 0.168, 0.198]
p_LTBI = [0.15, 0.2, 0.25]

F = lambda x: get_objective(x, ref, prm, gps_born, likelihood)



def MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, cov0, displ):
    d = len(x0)
    b = 0.05
    sd = sigma * (2.4**2) / d
    
    if fixinds is None:
        fixinds = []
        inds = fixinds[0]
        vals = fixinds[1]
    else:
        inds = []
        vals = []

    # Checks on the initial covariance matrix
    if cov0 is None:
        cov0 = np.eye(d)
        cov0[np.ix_(inds, inds)] = 0

    cov0[:blockind, blockind:] = 0
    cov0[blockind:, :blockind] = 0
    

    # Initialize output matrices
    xsto = np.zeros((d, n))
    outsto = np.zeros(n)
    history = np.zeros((d + 1, n))  # Rows: 1:d Proposed values 4. Accept or reject

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
            Y0 = multivariate_normal.rvs(X, 0.1**2 * cov0 * sigma / d)
            if t < 2 * d:
                Y = np.maximum(Y0, 0)
                Y[inds] = vals
            else:
                ind0 = 0
                ind1 = t - 1
                covmat = np.cov(xsto[:, ind0:ind1].T)
                covmat[np.ix_(inds, inds)] = 0
                covmat[:blockind, blockind:] = 0
                covmat[blockind:, :blockind] = 0
                covmat = (covmat + covmat.T) / 2

                # Precaution to make sure values don't get negative
                Y = np.maximum((1 - b) * multivariate_normal.rvs(X, sd * covmat) + b * Y0, 0)
                Y[inds] = vals
            history[0:d, t] = Y

            # --- Decide whether to accept or not
            FY = F(Y)
            if (npr.rand() < np.exp(FY - FX)) and (abs(FY) < np.inf):
                # Accept
                xsel = Y[:]
                FX = FY
                acc = acc + 1
                history[d, t] = 1
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

