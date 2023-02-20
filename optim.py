# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 22:13:20 2023

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.optimize import leastsq

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
    t7 = np.linspace(0,30000,30000+1)

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

    return (cInc[1:] - cInc[:-1])[-1]*100000 #this is the array for which i want the last value to be 7, and find the best beta for that, so this function returns the last value in that array
    

def residual(x):

    # Total population,  N.
    return np.sum((peak_infections(x) - 7) ** 2) #here i implemented what you said

x0 = 15 #init guess
res2 = minimize(residual, x0)
print(res2) 

 

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

    t7 = np.linspace(0,30000,30000+1) # KEY VARIATION

 

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

 

 

 

import pandas as pd

 

beta_inc = 0.5  # comment 16/0.5 = 32 increments

beta_range = [beta_inc * i for i in range(0, 40)]  # KEY VARIATION

cInc_range = [ round(peak_infections(beta_inc * i)*100000, 0) for i in range(0, 40)]

 

cIncidents = pd.DataFrame(np.column_stack([beta_range, cInc_range]), columns=['Beta values','Coincidents'])

 

print(cIncidents)















#############################################################################################

import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import math

import plotly.express as px

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

def optimise_beta(beta):
  N = 1 # Total population, N
  I0, R0 = 0.001, 0 # Initial number of infected and recovered individuals, I0 and R0.
  U0 = N - I0 - R0 # Everyone else, U0, is susceptible to infection initially.
  J0 = I0 
  Lf0, Ls0 = 0, 0

  # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
  gamma = 365/75
  int_gamma = 0.8
  mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
  u, v, w = 0.88, 0.083, 0.0006
  t = np.linspace(0, 500, 500+1)

  # Integrate the SIR equations over the time grid, t.
  solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
  U, Lf, Ls, I, R, cInc = solve.T
  J_diff = cInc[1:] - cInc[:-1]

  equilibrium_value = list(J_diff*100000)[-1]

  return equilibrium_value



beta_range = np.linspace(0,100,1000)
equilibrium_values = pd.DataFrame(columns=['beta','equilibrium_value'])


for idx, beta in enumerate(beta_range):
  equilibrium_values.at[idx,'beta'] = beta
  equilibrium_values.at[idx,'equilibrium_value'] = optimise_beta(beta)
  
  
equilibrium_values[(equilibrium_values['equilibrium_value']>=7) & (equilibrium_values['equilibrium_value']<=8)]  
  
fig = px.line(equilibrium_values, x="beta", y="equilibrium_value")
fig.show()

