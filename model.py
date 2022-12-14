# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:55:33 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
J0 = I0
Lf0, Ls0 = 0, 0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 12, 2/5
mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
u, v, w = 0.083, 0.88, 0.0006
t = np.linspace(0, 20, 20+1)

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
J_diff = np.diff(cInc)
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
    N = 100000
    # Initial number of infected and recovered individuals, I0 and R0.
    #I0, R0 = 10, 0
    # Everyone else, S0, is susceptible to infection initially.
    R0 = 0
    beta = x[0]
    I0 = x[1]
    gamma = x[2]
    S0 = N - I0 - R0
    J0 = I0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #reproductive no. R zero is beta/gamma
    #gamma = 1/6 #rate should be in weeks now
    # A grid of time points 
    times = np.arange(0,84,7)

    # The SIR model differential equations.
    def deriv(y, times, N, beta, gamma):
        S, I, R, J = y
        dS = ((-beta * S * I) / N)
        dI = ((beta * S * I) / N) - (gamma * I)
        dR = (gamma * I)
        dJ = ((beta * S * I) / N) #incidence
        return dS, dI, dR, dJ

    # Initial conditions are S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    solve = odeint(deriv, (S0, I0, R0, J0), times, args=(N, beta, gamma))
    S, I, R, J = solve.T

    return I/N

def residual(x, df):

    # Total population,  N.
    StartingPop = 100000
    incidence = df.incidence.to_numpy()/StartingPop
    return np.sum((peak_infections(x,df) - incidence) ** 2)


x0 = [0.5, 10, 1/6] #beta, i0, gamma
res = minimize(residual, x0, args=(df), method="Nelder-Mead", options={'fatol':1e-04}).x
print(res)

