# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 15:46:24 2023

@author: ISA
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta


def get_distribution_fns(data, distribution, show=True):
    prctiles = data  # Rename 'data' to 'prctiles' for clarity

    dat = sorted(prctiles)

    if distribution == 'lognorm':
        # Calculate parameters for the log-normal distribution
        mu = np.log(dat[1])
        sigma = (np.log(dat[2]) - np.log(dat[0])) / (2 * 1.96)  # Using 1.96 for the 95% interval
        params = (mu, sigma)

        # Define the log-pdf function for the log-normal distribution
        logfn = lambda x: -((np.log(x) - mu) ** 2) / (2 * sigma ** 2) - np.log(x * sigma * np.sqrt(2 * np.pi))

        # Adjust x values based on the mean
        x_mean = np.exp(mu)
        x = np.linspace(x_mean - 3 * sigma, x_mean + 3 * sigma, 100)

    elif distribution == 'beta':
        # Calculate parameters for the beta distribution
        a = 4 * dat[1] + 1  # Using the median for 'a'
        b = 4 * (1 - dat[1]) + 1  # Using (1 - median) for 'b'
        params = (a, b)

        # Define the log-pdf function for the beta distribution
        logfn = lambda x: (a - 1) * np.log(x) + (b - 1) * np.log(1 - x) - beta.logpdf(x, a, b)

        # Adjust x values based on the mean
        x_mean = a / (a + b)
        x = np.linspace(x_mean - 3 * np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))), 
                        x_mean + 3 * np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))), 100)

    else:
        raise ValueError('Invalid distribution type. Use "lognorm" or "beta".')

    if show:
        # Calculate the PDF values for the distribution
        pdf_values = np.exp([logfn(xi) for xi in x])

        # Create a plot for the PDF
        plt.figure(figsize=(8, 6))
        plt.plot(x, pdf_values, label=f'{distribution.capitalize()} PDF')
        plt.title(f'{distribution.capitalize()} Distribution PDF')
        plt.xlabel('X')
        plt.ylabel('PDF')
        plt.legend()
        plt.grid(True)
        plt.show()

    return
