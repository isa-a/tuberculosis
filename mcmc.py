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
from scipy.stats import norm  # Add this import
from scipy.stats import lognorm, beta

data = {
    'treatment_access': np.array([0.39, 0.41, 0.43]),
    'treatment_comp': np.array([0.76, 0.78, 0.8]),
    'mort': np.array([0.28, 0.3, 0.32]),
    'p_migrTB': np.array([0.708, 0.728, 0.748]),
    'p_migrpopn': np.array([0.138, 0.168, 0.198]),
    'p_LTBI': np.array([0.15, 0.2, 0.25])
} # define data lists



#takes data and fits it to distribution type
def get_distribution_fns(data, distribution_type, show=False):
    if distribution_type == 'lognorm':
        params = lognorm.fit(data, floc=0)
        dist = lognorm(*params)
    elif distribution_type == 'beta':
        params = beta.fit(data, floc=0, fscale=1)
        dist = beta(*params)
    else:
        raise ValueError("Unsupported distribution type")

    if show:
        # Create a range of x values for plotting
       x = np.linspace(data.min(), data.max(), 1000)

       # Plot the fitted distribution
       plt.plot(x, dist.pdf(x), label=distribution_type)
       plt.xlabel('Value')
       plt.ylabel('Probability Density')
       plt.legend()
       plt.title(f'{distribution_type} Distribution')

       # Show the plot
       plt.show()


    return dist

f1a = get_distribution_fns(data['treatment_access'], 'lognorm')
f1b = get_distribution_fns(data['treatment_comp'], 'lognorm')
f2 = get_distribution_fns(data['mort'], 'lognorm')
f3 = get_distribution_fns(data['p_migrTB'], 'beta')
f4 = get_distribution_fns(data['p_migrpopn'], 'beta')
f5 = get_distribution_fns(data['p_LTBI'], 'beta')


# f1a = get_distribution_fns(data['treatment_access'], 'lognorm', show=True)
# f1b = get_distribution_fns(data['treatment_comp'], 'lognorm', show=True)
# f2 = get_distribution_fns(data['mort'], 'lognorm', show=True)
# f3 = get_distribution_fns(data['p_migrTB'], 'beta', show=True)
# f4 = get_distribution_fns(data['p_migrpopn'], 'beta', show=True)
# f5 = get_distribution_fns(data['p_LTBI'], 'beta', show=True)

def likelihood(treatment_access, treatment_comp, mort, p_migrTB, p_migrpopn, p_LTBI):
    return f1a.pdf(treatment_access) + f1b.pdf(treatment_comp) + f2.pdf(mort) + f3.pdf(p_migrTB) + f4.pdf(p_migrpopn) + f5.pdf(p_LTBI)
# ... (previous code)

# MCMC settings
num_iterations = 100
beta_samples = []
gamma_samples = []

# Initialize current values
current_beta = np.random.normal(8, 2)
current_gamma = np.random.normal(0.4, 0.1)

# Define priors for beta and gamma
beta_prior_mean = 8
beta_prior_std = 2
gamma_prior_mean = 0.4
gamma_prior_std = 0.1


# MCMC loop
for i in range(num_iterations):
    # Calculate likelihood for current beta and gamma
    current_likelihood = likelihood(data['treatment_access'], data['treatment_comp'], data['mort'], data['p_migrTB'], data['p_migrpopn'], data['p_LTBI'])

    # Calculate log priors for beta and gamma
    current_log_beta_prior = norm.logpdf(current_beta, beta_prior_mean, beta_prior_std)
    current_log_gamma_prior = norm.logpdf(current_gamma, gamma_prior_mean, gamma_prior_std)

    # Calculate acceptance ratio
    acceptance_ratio = (current_likelihood + current_log_beta_prior + current_log_gamma_prior) - \
                       (np.log(np.random.rand()) + current_log_beta_prior + current_log_gamma_prior)

    # Accept or reject based on the acceptance ratio
    if np.any(acceptance_ratio > 0) or np.log(np.random.rand()) < acceptance_ratio:
        beta_samples.append(current_beta)
        gamma_samples.append(current_gamma)


# Calculate posterior means for beta and gamma
beta_mean = np.mean(beta_samples)
gamma_mean = np.mean(gamma_samples)

# Display the results
print("Estimated Beta:", beta_mean)
print("Estimated Gamma:", gamma_mean)








##################  ignore   ###########################
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

