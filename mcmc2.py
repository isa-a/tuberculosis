# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 21:21:49 2023

@author: ISA
"""

import numpy as np

def MCMC_adaptive(F, x0, n, sigma, cov0, displ=True):
    d = len(x0)
    b = 0.05
    sd = sigma * 2.4**2 / d

    # Checks on the initial covariance matrix
    if cov0 is None:
        cov0 = np.eye(d)

    xsto = np.zeros((n, d))
    outsto = np.zeros(n)
    history = np.zeros((n, d + 1))
    accept_rate = 0

    xsto[0] = x0
    xbar = xsto.copy()
    FX = F(x0)
    outsto[0] = FX
    acc = 0

    if displ:
        print("Starting MCMC...")

    for t in range(1, n):
        X = xsto[t - 1]

        # Make a proposal from the distribution
        Y0 = np.random.multivariate_normal(X, 0.1**2 * cov0 * sigma / d)
        Y = np.maximum(Y0, 0)

        history[t, :d] = Y

        FY = F(Y)
        if np.random.rand() < np.exp(FY - FX) and np.isfinite(FY):
            # Accept
            xsel = Y
            FX = FY
            acc += 1
            history[t, -1] = 1
        else:
            # Reject
            xsel = xsto[t - 1]

        xsto[t] = xsel
        outsto[t] = FX
        xbar[t] = (xbar[t - 1] * t + xsel) / (t + 1)

        # Display options
        if displ:
            print(f'Iteration {t} of {n} ({t/n*100:.1f}% complete)')

        
        
        # if displ and t % (n // 25) == 0:
        #     print(f"Progress: {t / n * 25:.2f}%")
        # if displ and t % 200 == 0:
        #     pass  # Plotting code here if needed

    accept_rate = acc / n

    return xsto, outsto, history, accept_rate



# Define your log-posterior density function F
def log_posterior(x):
    # Replace this with your actual log-posterior calculation
    return -0.5 * np.sum(x**2)

# Initial parameter set x0
x0 = [21.2257, 0.1499, 5.0908, 0.788, 21.493]

# Number of iterations
n = 1000

# Proposal standard deviation (sigma)
sigma = 0.1

# Optionally, specify the initial covariance matrix (cov0)
cov0 = np.eye(len(x0))

# Run the MCMC algorithm
xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, n, sigma, cov0)

# You can now analyze the results in the variables xsto, outsto, history, and accept_rate