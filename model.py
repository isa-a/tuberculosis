# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:55:33 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#three compartments, Susceptible S, infected I, recovered R
#dS/dt, dI/dt, dR/dt
#susceptible sees birth rate coming in, deaths leaving and force of infection leaving
#infected sees FOI coming in, deaths leaving and recovery rates
#recovered sees recovery rate coming in, deaths leaving
#beta is tranmission coefficient, FOI is beta * (I/N) where N is total pop
#initially consider a model not accounting for births and deaths


# Total population, N.
N = 100000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
Prevalence = I0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 1.47617188, 1/7
lamda = beta * I0
mu, muTB, sigma, rho = 1/80, 1/6, 1/6, 0.03
u, v, w = 0.083, 0.88, 0.0006

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    U, Lf, Ls, I, R, Prevalence = y
    b = mu(U + Lf + Ls + R) + muTB
    dU = b - ((lamda + mu) * U)
    dLf = (lamda*U) + ((c*lamda)*(Ls + R)) - ((u + v + mu) * Lf)
    dLs = (u * Lf) - ((w + c*lamda + mu) * Ls)
    dJ = ((beta * S * I) / N)
    return dS, dI, dR, dJ


# Total population, N.
N = 100000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
J0 = I0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 1.47617188, 1/7
# A grid of time points (in days)
t = np.linspace(0, 77, 77+1)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R, J = y
    dS = ((-beta * S * I) / N)
    dI = ((beta * S * I) / N) - (gamma * I)
    dR = (gamma * I)
    dJ = ((beta * S * I) / N)
    return dS, dI, dR, dJ

# Initial conditions are S0, I0, R0
# Integrate the SIR equations over the time grid, t.
solve = odeint(deriv, (S0, I0, R0, J0), t, args=(N, beta, gamma))
S, I, R, J = solve.T

J_diff = J[1:] - J[:-1]
J_diff = np.diff(J)
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
#ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2, label='Daily incidence')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.1)
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
plt.show()



#J_diff = J[1:] - J[:-1]
#J_diff = np.diff(J)
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S, 'b', alpha=1, lw=2, label='Susceptible')
#ax.plot(t, I, 'r', alpha=1, lw=2, label='Infected')
#ax.plot(t, R, 'black', alpha=1, lw=2, label='Recovered')
#ax.plot(t, J, 'green', alpha=1, lw=2, label='Incidence')
#ax.plot(t, J, 'red', alpha=1, lw=2, label='Cumulative incidence')
ax.plot(t[1:], J_diff, 'blue', alpha=1, lw=2, label='Daily incidence')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number')
#ax.set_ylim(0,1.1)
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
plt.show()