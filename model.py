# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:55:33 2022

@author: ISA
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import math
import pandas as pd


# Total population, N.
N = 1
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 0.001, 0
# Everyone else, S0, is susceptible to infection initially.
U0 = N - I0 - R0
J0 = I0
Lf0, Ls0 = 0, 0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta, gamma = 2,5
int_gamma = 0.8
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
#J_diffint = cIncint[1:] - cIncint[:-1]
#J_diff = np.diff(cInc)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, U*100000, 'black', alpha=1, lw=2, label='uninfected')
#ax.plot(t, Lf/100000, 'r', alpha=1, lw=2, label='latent fast')
#ax.plot(t, Ls/100000, 'black', alpha=1, lw=2, label='latent slow')
#ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
#ax.plot(t, R*100000, 'red', alpha=1, lw=2, label='recovered')
ax.plot(t[1:], J_diff*100000, 'blue', alpha=1, lw=2, label='incidence')
#ax.plot(t[1:]+2019, J_diffint*100000, 'red', alpha=1, lw=2, label='intervention incidence')
#ax.plot(t, cInc, 'red', alpha=1, lw=2, label='Prevalence')
ax.set_xlabel('Time in years')
ax.set_ylabel('Number')
#ax.set_xlim(2019, 2030)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#plt.title("Intervention")
plt.show()


# The SIR model differential equations.
def derivint(y, t, N, beta, int_gamma, mu, muTB, sigma, rho, u, v, w):
    U, Lf, Ls, I, R, cInc = y
    #int_gamma = t/2500+0.4
    b = (mu * (U + Lf + Ls + R)) + (muTB * I)
    lamda = beta * I
    clamda = 0.2 * lamda
    dU = b - ((lamda + mu) * U)
    dLf = (lamda*U) + ((clamda)*(Ls + R)) - ((u + v + mu) * Lf)
    dLs = (u * Lf) - ((w + clamda + mu) * Ls)
    dI = w*Ls + v*Lf - ((int_gamma + muTB + sigma) * I) + (rho * R)
    dR = ((int_gamma + sigma) * I) - ((rho + clamda + mu) * R)
    cI = w*Ls + v*Lf + (rho * R)
    return dU, dLf, dLs, dI, dR, cI


# Integrate the SIR equations over the time grid, t.
solveint = odeint(derivint, (U[-1], Lf[-1], Ls[-1], I[-1], R[-1], J0), t, args=(N, beta, int_gamma, mu, muTB, sigma, rho, u, v, w))
Uint, Lfint, Lsint, Iint, Rint, cIncint = solveint.T


def peak_infections(beta):
 
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

    def deriv(y, t7, N, beta, gamma, mu, muTB, sigma, rho, u, v, w):
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
    U, Lf, Ls, I, R, cInc = solve.T

    return (cInc[1:] - cInc[:-1])[-1]

def residual(x):

    # Total population,  N.
    return np.sum((peak_infections(x) - 7/100000) ** 2)

x0 = 9.5
res = minimize(residual, x0).x
print(res)


def error_function(params):
    gamma = params[0]
    beta = params[1]
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
    return (7/100000 - np.max(I)/N)**2

initial_guess = [0.55, 3]
result = minimize(error_function, initial_guess)
fit_params = result.x


J_diff = cInc[1:] - cInc[:-1]
#J_diffint = cIncint[1:] - cIncint[:-1]
#J_diff = np.diff(cInc)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, U*100000, 'black', alpha=1, lw=2, label='uninfected')
#ax.plot(t, Lf/100000, 'r', alpha=1, lw=2, label='latent fast')
#ax.plot(t, Ls/100000, 'black', alpha=1, lw=2, label='latent slow')
#ax.plot(t, I*100000, 'green', alpha=1, lw=2, label='infected')
#ax.plot(t, R*100000, 'red', alpha=1, lw=2, label='recovered')
ax.plot(t[1:], J_diff*100000, 'blue', alpha=1, lw=2, label='incidence')
#ax.plot(t[1:]+2019, J_diffint*100000, 'red', alpha=1, lw=2, label='intervention incidence')
#ax.plot(t, cInc, 'red', alpha=1, lw=2, label='Prevalence')
ax.set_xlabel('Time in years')
ax.set_ylabel('Number')
#ax.set_xlim(2019, 2030)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#plt.title("Intervention")
plt.show()


def lik(parameters):
    m = parameters[0]
    b = parameters[1]
    sigma = parameters[2]
    for i in np.arange(0, len(x)):
        y_exp = m * x + b
    L = (len(x)/2 * np.log(2 * np.pi) + len(x)/2 * np.log(sigma ** 2) + 1 /
         (2 * sigma ** 2) * sum((y - y_exp) ** 2))
    return L

def loglik(params):
    
    beta = params[0]
    gamma = params[1]
    
    solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
    U, Lf, Ls, I, R, cInc = solve.T #get trajectories

    muPrev, sigmaPrev = I[-1]*100000, 40 #I (prevalence)
    muInc, sigmaInc = (cInc[1:] - cInc[:-1])[-1]*100000, 30 #cInc (incidence)
    #n = 10000
    
    # logPrev = np.random.lognormal(np.log((muPrev**2) / (muPrev**2 + sigmaPrev**2)**0.5), (np.log(1 + (sigmaPrev**2 / muPrev**2)))**0.5, n) #lognormal
    # logInc = np.random.lognormal(np.log((muInc**2) / (muInc**2 + sigmaInc**2)**0.5), (np.log(1 + (sigmaInc**2 / muInc**2)))**0.5, n) #lognormal
    
    xPrev = I[-1]*100000 #value of x in formula for log of pdf
    xInc = (cInc[1:] - cInc[:-1])[-1]*100000 #value of x in formula for log of pdf
    
    logmuPrev = np.log((muPrev**2) / (muPrev**2 + sigmaPrev**2)**0.5) #lognormal params
    logsdPrev = (np.log(1 + (sigmaPrev**2 / muPrev**2)))**0.5
    
    logmuInc = np.log((muInc**2) / (muInc**2 + sigmaInc**2)**0.5)#lognormal params
    logsdInc = (np.log(1 + (sigmaInc**2 / muInc**2)))**0.5
    
    L_prev = -0.5*((np.log(xPrev) - logmuPrev) / logsdPrev)**2 - np.log(xPrev * logsdPrev * (2*math.pi)**0.5) #log of pdf for prev and inc
    L_inc = -0.5*((np.log(xInc) - logmuInc) / logsdInc)**2 - np.log(xInc * logsdInc * (2*math.pi)**0.5)
    
    logsum = L_prev + L_inc #summing logs
    np.exp(logsum) #exp for likelihood
    return np.exp(logsum)


x = I[:-1]
y = cInc[1:] - cInc[:-1]
x0 = [8, 0.4]
lik_model = minimize(loglik, x0)

mu, sigma = 5.69, 0.01 # mean and standard deviations 
s= np.random.lognormal(mu, sigma, 1000)
import matplotlib.pyplot as plt
count, bins, ignored = plt.hist(s, 100, density=True, align='mid')
x = np.linspace(min(bins), max(bins), 10000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
       / (x * sigma * np.sqrt(2 * np.pi)))
plt.plot(x, pdf, linewidth=2, color='r')
plt.axis('tight')
plt.show()


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


