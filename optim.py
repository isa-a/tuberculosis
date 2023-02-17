# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 22:13:20 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import scipy.optimize

def peak_infections(beta): #contains initial values, but need to estimate beta
 
    N = 1
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 0.001, 0
    # Everyone else, S0, is susceptible to infection initially.
    U0 = N - I0 - R0
    J0 = I0
    Lf0, Ls0 = 0, 0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    gamma = 365/75
    mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
    u, v, w = 0.88, 0.083, 0.0006
    t7 = np.linspace(0,500,500+1)

    def deriv(y, t7, N, beta, gamma, mu, muTB, sigma, rho, u, v, w): #solves ode system
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
    solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t7, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
    U, Lf, Ls, I, R, cInc = solve.T #output trajectories

    return (cInc[1:] - cInc[:-1])[-1] #this is the array for which i want the last value to be 7, and find the best beta for that, so this function returns the last value in that array

def residual(x):

    # Total population,  N.
    return np.sum((peak_infections(x) - 7/100000) ** 2) #here i implemented what you said

x0 = 9.5 #init guess
res = minimize
print(res) #just returns initial guess :(

