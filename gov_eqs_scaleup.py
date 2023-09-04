# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:33:40 2023

@author: ISA
"""

from model_setup import *
from mcmc import beta_mean, gamma_mean
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the start and end years
start = 2000
end = 2020

# Calculate the number of time steps based on the desired range
num_years = end_year - start_year + 1  # +1 to include the end year
num_time_steps = num_years * 365  # Assuming daily time steps

# Create the time array to represent years from 2010 to 2020
time_range = np.linspace(start, end, num_time_steps)

solve = odeint(gov_eqs, (U0, Lf0, Ls0, I0, R0), time_range, args=(N, beta_mean, gamma_mean, u, v, w))


fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, U*100000, 'black', alpha=1, lw=2, label='uninfected')
#ax.plot(t, Lf*100000, 'black', alpha=1, lw=2, label='latent fast')
#ax.plot(t, Ls*100000, 'purple', alpha=1, lw=2, label='latent slow')
ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
ax.plot(t, R*100000, 'blue', alpha=1, lw=2, label='recovered')
ax.set_xlabel('Time')
ax.set_ylabel('Number')
#ax.set_xlim(2019, 2030)
ax.grid(which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#plt.title("Incidence")
plt.show()
