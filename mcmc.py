# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 16:58:31 2023

@author: ISA
"""

#MCMC to infer params from model_setup


from model_setup import *
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from tqdm import tqdm

# Define the observed data
observed_data = I  # uses observed data for the number for I

# Define the prior distributions for beta and gamma
beta_prior_mean = 8  # mean beta
beta_prior_std = 2  # sd beta
gamma_prior_mean = 0.4  #  mean  gamma
gamma_prior_std = 0.1  #  sd gamma


#sds for the proposal distributions of beta and gamma
#dist generates new values for beta and gamma during mcmc
beta_proposal_std = 0.5  # Standard deviation for the proposal distribution of beta
gamma_proposal_std = 0.05  # Standard deviation for the proposal distribution of gamma

#num of iterations for mcmc
num_iterations = 10

beta_lower_bound = 0.0
beta_upper_bound = 10.0
gamma_lower_bound = 0.0
gamma_upper_bound = 1.0

#lists to store param samples
beta_samples = []
gamma_samples = []
#initialize the current values of beta and gamma as random samples from 
#the dists
current_beta = np.random.normal(beta_prior_mean, beta_prior_std)
current_gamma = np.random.normal(gamma_prior_mean, gamma_prior_std)

# calc likelihood and prior for the current beta and gamma
solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), t, args=(N, current_beta, current_gamma, u, v, w))
I_current = solve.T[3]  #get  infected individuals from the model 
current_likelihood = np.exp(-0.5 * np.sum((I_current - observed_data) ** 2))
current_beta_prior = np.exp(-0.5 * ((current_beta - beta_prior_mean) / beta_prior_std) ** 2)
current_gamma_prior = np.exp(-0.5 * ((current_gamma - gamma_prior_mean) / gamma_prior_std) ** 2)


progress_bar = tqdm(total=num_iterations, ncols=80, bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')

# perform the mcmc
#in the loop new values for params are proposed from distributions using normal dist
for i in range(num_iterations): #iterates num_iterations
    proposed_beta = np.random.normal(current_beta, beta_proposal_std)
    proposed_gamma = np.random.normal(current_gamma, gamma_proposal_std)
    
# Check if the proposed beta and gamma values are within the specified bounds
    if beta_lower_bound <= proposed_beta <= beta_upper_bound and \
            gamma_lower_bound <= proposed_gamma <= gamma_upper_bound:
        # Update the model with the proposed beta and gamma values
        solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), t, args=(N, proposed_beta, proposed_gamma, u, v, w))
        I_proposed = solve.T[3]  # Get the infected individuals from the model output

        # Calculate the likelihood of the proposed beta and gamma values
        proposed_likelihood = np.exp(-0.5 * np.sum((I_proposed - observed_data) ** 2))

        # Calculate the prior probabilities of the proposed beta and gamma values
        proposed_beta_prior = np.exp(-0.5 * ((proposed_beta - beta_prior_mean) / beta_prior_std) ** 2)
        proposed_gamma_prior = np.exp(-0.5 * ((proposed_gamma - gamma_prior_mean) / gamma_prior_std) ** 2)

        # Calculate the acceptance ratio
        acceptance_ratio = (proposed_likelihood * proposed_beta_prior * proposed_gamma_prior) / \
                           (current_likelihood * current_beta_prior * current_gamma_prior + 1e-10)    
    # Accept or reject the proposed beta and gamma values based on the acceptance ratio
    if np.random.rand() < acceptance_ratio:
        current_beta = proposed_beta
        current_gamma = proposed_gamma
        current_likelihood = proposed_likelihood
        current_beta_prior = proposed_beta_prior
        current_gamma_prior = proposed_gamma_prior
    
    # Store the samples
    beta_samples.append(current_beta)
    gamma_samples.append(current_gamma)
    
    # Update the progress bar
    progress_bar.update(1)

# Close the progress bar
progress_bar.close()

# Compute the mean values of beta and gamma
beta_mean = np.mean(beta_samples)
gamma_mean = np.mean(gamma_samples)

# Print the results
print("Inferred values:")
print("Beta:", beta_mean)
print("Gamma:", gamma_mean)

# Plot the posterior distributions of beta and gamma
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.hist(beta_samples, bins=30, density=True, color='blue', alpha=0.7)
plt.axvline(x=beta_mean, color='red', linestyle='--', label='Mean')
plt.xlabel('Beta')
plt.ylabel('Posterior Density')
plt.title('Posterior Distribution of Beta')
plt.legend()
##########################################################################
plt.subplot(2, 1, 2)
plt.hist(gamma_samples, bins=30, density=True, color='green', alpha=0.7)
plt.axvline(x=gamma_mean, color='red', linestyle='--', label='Mean')
plt.xlabel('Gamma')
plt.ylabel('Posterior Density')
plt.title('Posterior Distribution of Gamma')
plt.legend()

plt.tight_layout()
plt.show()


# traces
plt.figure(figsize=(8, 4))
plt.plot(range(num_iterations), beta_samples, color='blue', alpha=0.7)
plt.axhline(y=beta_mean, color='red', linestyle='--', label='Mean')
plt.xlabel('Iteration')
plt.ylabel('Beta')
plt.title('Trace Plot for Beta')
plt.legend()
plt.show()
########################################################################
plt.figure(figsize=(8, 4))
plt.plot(range(num_iterations), gamma_samples, color='green', alpha=0.7)
plt.axhline(y=gamma_mean, color='red', linestyle='--', label='Mean')
plt.xlabel('Iteration')
plt.ylabel('Gamma')
plt.title('Trace Plot for Gamma')
plt.legend()
plt.show()