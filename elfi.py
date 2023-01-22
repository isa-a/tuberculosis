# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 14:16:08 2023

@author: ISA
"""
import sys
import elfi
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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
beta0, gamma0 = 8, 0.4
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

def derivative(beta, gamma, batch_size = 1, random_state = None):
    
    y0 = [U0, Lf0, Ls0, I0, R0, J0]
    
    times = np.linspace(0, 500, 500+1)
    
    resultz = odeint(deriv, y0, times, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))

    return resultz

# def SSE(x, y):
#     return np.atleast_1d(np.sum((x-y) ** 2, axis=1))

def select_var(X, variable=3):
    """
    variable name
       0       U0
       1      Lf0
       2      Ls0
       3       I0
       4       R0
       5       J0
    """
    return X[:,:,variable]

# def autocov(x, lag=1):
#     x = np.atleast_2d(x)
#     # In R this is normalized with x.shape[1]
#     C = np.mean(x[:, lag:] * x[:, :-lag], axis=1)
#     return C

vectorized_derivative = elfi.tools.vectorize(derivative)

y_obs_I = vectorized_derivative(beta0, gamma0)

model = elfi.new_model()

beta_prior = elfi.Prior('uniform', 7, 9, model=model)
gamma_prior = elfi.Prior('uniform', 0, 1, model=model)

sim_results = elfi.Simulator(vectorized_derivative, model['beta_prior'], model['gamma_prior'], observed = y_obs_I)
elfi.draw(sim_results)
sim_results.generate()

var = elfi.Summary(select_var, model['sim_results'], 3)

#summ = elfi.Summary(autocov, model['I'], 1)

#elfi.Distance(SSE, model['S1'], name='dist')
elfi.Distance('euclidean', var, name='dista')
d = elfi.Distance('euclidean', var)


rej = elfi.Rejection(model['dista'], batch_size = 1, seed = 1)
result = rej.sample(1000, quantile = 0.01)


fig, ax = plt.subplots()
rej.plot_state(ax = ax)
ax.set_xlim([0, 15])
ax.set_ylim([0, 2])
plt.savefig("parameter_space.png")
plt.close()


fig, ax = plt.subplots()
rej.plot_traces(ax = ax)
ax.set_xlim([0, 15])
ax.set_ylim([0, 2])
plt.savefig("tra.png")
plt.close()


rej.extract_result()

bolfi = elfi.BOLFI(d, batch_size=1, initial_evidence=100, update_interval=10, bounds={'beta_prior':(6, 9), 'gamma_prior':(0, 1)}, seed=1)
post = bolfi.fit(n_evidence=500)


bolfi.plot_discrepancy()
plt.savefig("bolfi_discrepancy.pdf")
plt.close()

sys.stderr.write("Sampling from BOLFI posterior\n")
result_BOLFI = bolfi.sample(10000, algorithm="metropolis")

result_BOLFI.plot_traces()
plt.savefig("posterior_traces.pdf")
plt.close()

post2 = bolfi.extract_posterior(-1.)
post.plot(logpdf=True)
