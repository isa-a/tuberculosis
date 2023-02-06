# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:03:57 2023

@author: ISA
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math
import scipy.stats as sps
from scipy.integrate import odeint
import model

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

solve = odeint(model.deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
U, Lf, Ls, I, R, cInc = solve.T
#print(cInc[1:] - cInc[:-1])

muPrev, sigmaPrev = 400, 40 #I
muInc, sigmaInc = 300, 30 #cInc
n = 10000

logPrev = np.random.lognormal(np.log((muPrev**2) / (muPrev**2 + sigmaPrev**2)**0.5), (np.log(1 + (sigmaPrev**2 / muPrev**2)))**0.5, n) #lognormal
logInc = np.random.lognormal(np.log((muInc**2) / (muInc**2 + sigmaInc**2)**0.5), (np.log(1 + (sigmaInc**2 / muInc**2)))**0.5, n) #lognormal

xPrev = I[-1]*100000
xInc = (cInc[1:] - cInc[:-1])[-1]*100000

logmuPrev = np.log((muPrev**2) / (muPrev**2 + sigmaPrev**2)**0.5)
logsdPrev = (np.log(1 + (sigmaPrev**2 / muPrev**2)))**0.5

logmuInc = np.log((muInc**2) / (muInc**2 + sigmaInc**2)**0.5)
logsdInc = (np.log(1 + (sigmaInc**2 / muInc**2)))**0.5

L_prev = -0.5*((np.log(xPrev) - logmuPrev) / logsdPrev)**2 - np.log(xPrev * logsdPrev * (2*math.pi)**0.5)
L_inc = -0.5*((np.log(xInc) - logmuInc) / logsdInc)**2 - np.log(xInc * logsdInc * (2*math.pi)**0.5)

logsum = L_prev + L_inc
np.exp(logsum)







muPrev, sigmaPrev = 400, 40.
a = np.random.normal(muPrev, sigmaPrev, 100000)
count, bins, ignored = plt.hist(a, 1000, density=True)
plt.plot(bins, 1/(sigmaPrev * np.sqrt(2 * np.pi)) *
           np.exp( - (bins - muPrev)**2 / (2 * sigmaPrev**2) ),linewidth=3, color='r')
x = np.percentile(a, 2.5) # return 50th percentile, e.g median.
y = np.percentile(a, 97.5) # return 50th percentile, e.g median.
z = np.percentile(a, 50) # return 50th percentile, e.g median.
print(x)
print(y)
log_norm = np.exp(a)





muInc, sigmaInc = 300, 30
s = np.random.normal(muInc, sigmaInc, 100000)
count, bins, ignored = plt.hist(s, 1000, density=True)
plt.plot(bins, 1/(sigmaInc * np.sqrt(2 * np.pi)) *
           np.exp( - (bins - muInc)**2 / (2 * sigmaInc**2) ),linewidth=3, color='r')
p = np.percentile(s, 2.5) # return 50th percentile, e.g median.p = np.percentile(s, 2.5) # return 50th percentile, e.g median.
q = np.percentile(s, 97.5) # return 50th percentile, e.g median.
print(p)
print(q)


mu, sd = 300, 30
n = 100_000

# draw samples from distributions
a = np.random.normal(mu, sd, n)
b = np.random.lognormal(np.log((mu**2) / (mu**2 + sd**2)**0.5), (np.log(1 + (sd**2 / mu**2)))**0.5, n)

# use Scipy for analytical PDFs
d1 = sps.norm(mu, sd)
# warning: scipy parameterises its distributions very strangely
d2 = sps.lognorm(sd / mu, scale=mu)

# bins to use for histogram and x for PDFs
lo, hi = np.min([a, b]), np.max([a, b])
dx = (hi - lo) * 0.06
bins = np.linspace(lo, hi, 101)
x = np.linspace(lo - dx, hi + dx, 501)

# draw figure
fig, [ax1, ax2] = plt.subplots(nrows=2, sharex=True, sharey=True, figsize=(8, 5))

ax1.set_title("Incidence Normal draws")
ax1.set_xlim(lo - dx, hi + dx)
ax1.hist(a, bins, density=True, alpha=0.5)
ax1.plot(x, d1.pdf(x))
ax1.plot(x, d2.pdf(x), '--')

ax2.set_title("Log-Normal draws")
ax2.hist(b, bins, density=True, alpha=0.5)
ax2.plot(x, d1.pdf(x), '--', label="Normal PDF")
ax2.plot(x, d2.pdf(x), label="Log-Normal PDF")

ax2.legend()
fig.supylabel("Density")

x=1
mu=5.698807309229617
sd=0.0997513451195927
func = -0.5*((np.log(300) - mu) / sd)**2 - np.log(300 * sd * (2*math.pi)**0.5)

scipy.stats.norm.pdf(6,2.0,1.0)
print( np.log(scipy.stats.norm.pdf(a,2.5,97.5)).sum() )

x = np.linspace(-600, 600, 1000, endpoint=True)
y = []
for i in x:
    y.append(np.log(scipy.stats.norm.pdf(a,i,40)).sum())

plt.plot(x,y)
plt.title(r'Log-Likelihood')
plt.xlabel(r'$\mu$')

plt.grid()

print('mean ---> ', np.mean(a))
y_min = y.index(max(y))
print('mean (from max log likelohood) ---> ', x[y_min])


import numpy as np
import math
import scipy.optimize as optimize

def llnorm(mu, sigma, data):
    #mu, sigma = par
    ll = np.sum(math.log(2*math.pi*(sigma**2))/2 + ((data-mu)**2)/(2 * (sigma**2)))
    return ll
data = a
result = optimize.minimize(llnorm, [2,1], args = (data))
mu, sigma = 300, 30
llnorm(300, 30, s)
math.log(482098)

mu, sigma = 400, 40. # mean and standard deviation
q = np.random.lognormal(mu, sigma, 1000)
count, bins, ignored = plt.hist(q, 100, density=True, align='mid')
x = np.linspace(min(bins), max(bins), 10000)
pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
       / (x * sigma * np.sqrt(2 * np.pi)))
plt.plot(x, pdf, linewidth=2, color='r')
plt.axis('tight')
plt.show()
