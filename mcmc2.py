# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 21:21:49 2023

@author: ISA
"""

import numpy as np
from obj import get_objective
from setup_model import r,p,agg,sel,ref,xi,prm,gps_born,likelihood_function
from tqdm import tqdm
import time

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
    #accept_rate = 0

    xsto[0] = x0
    xbar = xsto.copy()
    FX = F(x0)
    outsto[0] = FX
    acc = 0

    if displ:
        print("Starting MCMC...")

    start_time = time.time()  # Start the timer

    #for t in range(1, n):
    # Wrapping the range with tqdm for progress bar display
    for t in tqdm(range(1, n), total=n-1, desc="MCMC Progress", ncols=100):
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
        
        end_time = time.time()  # End the timer

        elapsed_time = end_time - start_time  # Calculate elapsed time


        # Display options
        if displ:
            #print('Iteration', t, 'of', n, '({:.1f}% complete)'.format(t/n*100))
            print(f"MCMC completed in {elapsed_time:.2f} seconds.")

    accept_rate = acc / n

    return xsto, outsto, history, accept_rate



# obj = lambda x: get_objective(x, ref, prm, gps_born,likelihood_function)[0]

# # Initial parameter set x0
# x0 = [21.2257, 0.1499, 5.0908, 0.788, 21.493]

# # Number of iterations
# n = 50

# # Proposal standard deviation (sigma)
# sigma = 0.1

# # Optionally, specify the initial covariance matrix (cov0)
# cov0 = np.eye(len(x0))

# # Initial MCMC run
# xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, n, sigma, cov0)

# # Find the parameter set with the maximum log-posterior density
# inds = np.where(outsto == np.max(outsto))[0]
# x0 = xsto[inds[0], :]

# # Run MCMC again with updated cov0 (without blockinds or fixinds)
# cov0 = np.cov(xsto.T)
# xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, n, sigma, cov0)












































# Run the MCMC algorithm
# xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x0, n, sigma, cov0)














# import matplotlib.pyplot as plt

# # Assuming 'samples' is a NumPy array with shape (n_iterations, n_parameters)
# n_iterations, n_parameters = xsto.shape

# # Create a figure with subplots for each parameter
# fig, axes = plt.subplots(n_parameters, figsize=(8, 2 * n_parameters), sharex=True, dpi=800)

# # Loop over each parameter and plot its trace
# for i in range(n_parameters):
#     ax = axes[i]
#     parameter_trace = xsto[:, i]  # Extract the trace for the i-th parameter
#     ax.plot(parameter_trace, color='b', alpha=0.7, label=f'Parameter {i + 1}')
#     ax.set_ylabel(f'Parameter {i + 1}')
#     ax.axhline(y=parameter_trace.mean(), color='r', linestyle='--', label='Mean')
#     ax.legend()

# # Add a common x-axis label (assuming the x-axis represents iterations)
# axes[-1].set_xlabel('Iteration')

# # Adjust spacing between subplots
# plt.tight_layout()

# # Show the trace plot
# plt.show()

