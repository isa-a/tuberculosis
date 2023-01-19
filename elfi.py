# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 12:29:50 2023

@author: ISA
"""

import elfi
import time
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import logging
from scipy.integrate import odeint
from scipy.optimize import minimize
elfi.set_client('multiprocessing')

seed = 20170530  # this will be separately given to ELFI
np.random.seed(seed)

# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
J0 = I0
Lf0, Ls0 = 0, 0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 8, 0.4
mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
u, v, w = 0.88, 0.083, 0.0006
t = np.linspace(0, 500, 500+1)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma, mu, muTB, sigma, rho, u, v, w):
    U, Lf, Ls, I, R, cInc = y
    b = (mu * (U + Lf + Ls + R)) + (muTB * I)
    lamda = beta * I
    clamda = 0.2 * lamda
    dU = b - ((lamda + mu) * U)
    dLf = (lamda*U) + ((clamda)*(Ls + R)) - ((u + v + mu) * Lf)
    dLs = (u * Lf) - ((w + clamda + mu) * Ls)
    dI = w*Ls + v*Lf - ((gamma + muTB + sigma) * I) + (rho * R)
    dR = ((gamma + sigma) * I) - ((rho + clamda + mu) * R)
    cI = w*Ls + v*Lf + (rho * R)
    return dU, dLf, dLs, dI, dR, cI


# Integrate the SIR equations over the time grid, t.
solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
U, Lf, Ls, I, R, cInc = solve.T


#try using I as observed data
y_obs_I = np.array(I)
y_obs_I = y_obs_I.flatten()


def derivative(beta, gamma, batch_size = 1, random_state = None):
    
    y0 = [U0, Lf0, Ls0, I0, R0, J0]
    
    times = np.linspace(0, 500, 500+1)
    
    resultz = odeint(deriv, y0, times, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
    
    return resultz.flatten()


beta_prior = elfi.Prior('uniform', 6, 9, name = 'beta')
gamma_prior = elfi.Prior('uniform', 0, 1, name = 'gamma')

Y_lv = elfi.Simulator(derivative, beta_prior, gamma_prior, observed = y_obs_I)
Y_lv.generate()

# Define sum of squared errors as the distance function
def SSE(x, y):
   return np.atleast_1d(np.sum((x-y)**2))
d_lv = elfi.Distance(SSE, Y_lv)
d_lv.generate()
elfi.draw(d_lv)
rej = elfi.Rejection(d_lv, batch_size = 1, seed = seed)
result = rej.sample(50, quantile = 0.01)


fig, ax = plt.subplots()
rej.plot_state(ax = ax)
ax.set_xlim([0, 11])
ax.set_ylim([0, 1])
plt.savefig("parameter_space.png")
plt.close()



rej.extract_result()


result.summary()
import matplotlib.pyplot as plt
result.plot_pairs();
plt.show()
