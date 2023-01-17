# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:55:33 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize


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

J_diff = cInc[1:] - cInc[:-1]
#J_diff = np.diff(cInc)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, U/100000, 'b', alpha=1, lw=2, label='uninfected')
#ax.plot(t, Lf/100000, 'r', alpha=1, lw=2, label='latent fast')
#ax.plot(t, Ls/100000, 'black', alpha=1, lw=2, label='latent slow')
ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
#ax.plot(t, R/100000, 'red', alpha=1, lw=2, label='recovered')
ax.plot(t[1:], J_diff*100000, 'blue', alpha=1, lw=2, label='Daily incidence')
#ax.plot(t, cInc, 'red', alpha=1, lw=2, label='Prevalence')
ax.set_xlabel('Time in years')
ax.set_ylabel('Number')
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()


beta_samples = np.random.uniform(0, 30, 10)
gamma_samples = np.random.uniform(0, 2, 10)

for i in beta_samples:
    # Total population, N.
    N = 1
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 0.001, 0
    # Everyone else, S0, is susceptible to infection initially.
    U0 = N - I0 - R0
    J0 = I0
    Lf0, Ls0 = 0, 0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    beta, gamma = i, 0.4
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
    
    rejected = []
    accepted = []
    
    if 320 < I[-1]*100000 < 480 and 240 < (cInc[1:] - cInc[:-1])[-1]*100000 < 360:
        print('pprevalence is ', I[-1]*100000, 'incidence is ', (cInc[1:] - cInc[:-1])[-1]*100000)
    else:
        print('values rejected')

    
    J_diff = cInc[1:] - cInc[:-1]
    #J_diff = np.diff(cInc)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    #ax.plot(t, U/100000, 'b', alpha=1, lw=2, label='uninfected')
    #ax.plot(t, Lf/100000, 'r', alpha=1, lw=2, label='latent fast')
    #ax.plot(t, Ls/100000, 'black', alpha=1, lw=2, label='latent slow')
    ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
    #ax.plot(t, R/100000, 'red', alpha=1, lw=2, label='recovered')
    ax.plot(t[1:], J_diff*100000, 'blue', alpha=1, lw=2, label='Daily incidence')
    #ax.plot(t, cInc, 'red', alpha=1, lw=2, label='Prevalence')
    ax.set_xlabel('Time in years')
    ax.set_ylabel('Number')
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    plt.show()





t = np.linspace(0, 77, 77+1)

def peak_infections(x):
 
    # Weeks for which the ODE system will be solved
    #weeks = df.Week.to_numpy()

    # Total population, N.
    N = 1
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 0.001, 0
    # Everyone else, S0, is susceptible to infection initially.
    beta = x[0]
    gamma = x[1]
    U0 = N - I0 - R0
    J0 = I0
    Lf0, Ls0 = 0, 0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    beta, gamma = 15, 2/5
    mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
    u, v, w = 0.083, 0.88, 0.0006

    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    #gamma = 1/6 #rate should be in weeks now
    # A grid of time points 
    times = np.arange(0,20,2.5)

    # The SIR model differential equations.
    def deriv(y, times, N, beta, gamma, mu, muTB, sigma, rho, u, v, w):
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

    # Initial conditions are S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), times, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
    U, Lf, Ls, I, R, cInc = solve.T

    return I/N

def residual(x):

    # Total population,  N.
    StartingPop = 1
    prev = 0.003/StartingPop
    return np.sum((peak_infections(x) - prev) ** 2)


x0 = [6, 0.4] #beta, i0, gamma
res = minimize(residual, x0, method="Nelder-Mead", options={'fatol':1e-04}).x
print(res)

