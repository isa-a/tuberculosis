# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 19:11:28 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize

beta_samples = np.random.uniform(0, 30, 500)
gamma_samples = np.random.uniform(0, 2, 500)
accepted = []
rejected = []
for i, j in zip(beta_samples, gamma_samples):
    # Total population, N.
    N = 1
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 0.001, 0
    # Everyone else, S0, is susceptible to infection initially.
    U0 = N - I0 - R0
    J0 = I0
    Lf0, Ls0 = 0, 0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    beta, gamma = i, j
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

    if 320 < I[-1]*100000 < 480 and 240 < (cInc[1:] - cInc[:-1])[-1]*100000 < 360:
        acc = [320 < I[-1]*100000 < 480]
        accepted.append(acc)
        print('for beta of', beta, 'and gamma of', gamma, 'pprevalence is ', I[-1]*100000, 'incidence is ', (cInc[1:] - cInc[:-1])[-1]*100000)
    else:
        rejected.append(beta_samples)
        print('values of', beta, 'and gamma of', gamma, 'rejected')
print(len(accepted), 'vals accepted')
