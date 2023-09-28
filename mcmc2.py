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





import matplotlib.pyplot as plt

# Assuming 'samples' is a NumPy array with shape (n_iterations, n_parameters)
n_iterations, n_parameters = xsto.shape

# Create a figure with subplots for each parameter
fig, axes = plt.subplots(n_parameters, figsize=(8, 2 * n_parameters), sharex=True, dpi=800)

# Loop over each parameter and plot its trace
for i in range(n_parameters):
    ax = axes[i]
    parameter_trace = xsto[:, i]  # Extract the trace for the i-th parameter
    ax.plot(parameter_trace, color='b', alpha=0.7, label=f'Parameter {i + 1}')
    ax.set_ylabel(f'Parameter {i + 1}')
    ax.axhline(y=parameter_trace.mean(), color='r', linestyle='--', label='Mean')
    ax.legend()

# Add a common x-axis label (assuming the x-axis represents iterations)
axes[-1].set_xlabel('Iteration')

# Adjust spacing between subplots
plt.tight_layout()

# Show the trace plot
plt.show()

