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
from scipy.stats import norm

# defining the observed data
observed_data = I  # uses observed data for the number for I

# defining stats for priors of params
beta_prior_mean = 8  # mean beta
beta_prior_std = 2  # sd beta
gamma_prior_mean = 0.4  #  mean  gamma
gamma_prior_std = 0.1  #  sd gamma


#sds for the proposal distributions of beta and gamma
#during each iteration of the MCMC new values for 
#params are proposed by sampling from these dists
#proposed values are evaluated based on their likelihood and priors
beta_proposal_std = 0.75  # Standard deviation for the proposal distribution of beta
gamma_proposal_std = 0.2  # Standard deviation for the proposal distribution of gamma

#num of iterations for mcmc
num_iterations = 100

#param bounds
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

# calc log likelihood and prior for the current beta and gamma
#model equations are solved using the current values of params
#the resulting infected individuals are compared to the observed data
#to calculate log likelihood. log priors probabilities of the current params values are also computed
solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), t, args=(N, current_beta, current_gamma, u, v, w))
I_current = solve.T[3]  # Get infected individuals from the model
#calculate current LL 
current_log_likelihood = -0.5 * np.sum((I_current - observed_data) ** 2) - 0.5 * len(observed_data) * np.log(2 * np.pi)

# Use log-normal pdf to calculate log priors of current params
current_log_beta_prior = -0.5 * ((current_beta - beta_prior_mean) / beta_prior_std) ** 2 - np.log(beta_prior_std * np.sqrt(2 * np.pi))
current_log_gamma_prior = -0.5 * ((current_gamma - gamma_prior_mean) / gamma_prior_std) ** 2 - np.log(gamma_prior_std * np.sqrt(2 * np.pi))



progress_bar = tqdm(total=num_iterations, ncols=80, bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')

# Perform the MCMC
# Loop iterates num of samples
for i in range(num_iterations):
    # new vals for params are proposed from their distributions
    proposed_beta = np.random.normal(current_beta, beta_proposal_std)
    proposed_gamma = np.random.normal(current_gamma, gamma_proposal_std)

    # Model is updated with the new proposed vals and the log likelihood and log priors are calculated
    solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), t, args=(N, proposed_beta, proposed_gamma, u, v, w))
    I_proposed = solve.T[3]

    proposed_log_likelihood = -0.5 * np.sum((I_proposed - observed_data) ** 2) - 0.5 * len(observed_data) * np.log(2 * np.pi)
    proposed_log_beta_prior = -0.5 * ((proposed_beta - beta_prior_mean) / beta_prior_std) ** 2 - np.log(beta_prior_std * np.sqrt(2 * np.pi))
    proposed_log_gamma_prior = -0.5 * ((proposed_gamma - gamma_prior_mean) / gamma_prior_std) ** 2 - np.log(gamma_prior_std * np.sqrt(2 * np.pi))

    #acceptance ratio is calculated based on the proposed and current values of
    # params, log likelihood, and log prior probabilities
    acceptance_ratio = (proposed_log_likelihood + proposed_log_beta_prior + proposed_log_gamma_prior) - \
                       (current_log_likelihood + current_log_beta_prior + current_log_gamma_prior)

    # proposed vals are evaluated based on the acceptance ratio
    # and current vals are updated
    if np.log(np.random.rand()) < acceptance_ratio:
        current_beta = proposed_beta
        current_gamma = proposed_gamma
        current_log_likelihood = proposed_log_likelihood
        current_log_beta_prior = proposed_log_beta_prior
        current_log_gamma_prior = proposed_log_gamma_prior

    # Current vals stored in lists
    beta_samples.append(current_beta)
    gamma_samples.append(current_gamma)
    progress_bar.update(1)

# get mean values of params
beta_mean = np.mean(beta_samples)
gamma_mean = np.mean(gamma_samples)

parameter_list = [beta_mean, gamma_mean]

print("Inferred values:")
print("Beta:", beta_mean)
print("Gamma:", gamma_mean)


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

