# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:03:57 2023

@author: ISA
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math



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


#dist = np.sqrt( (x - (2.5/100))**2 + (y - (97.5/100))**2 + (z - (50/100))**2)



muInc, sigmaInc = 300, 30
s = np.random.normal(muInc, sigmaInc, 100000)
count, bins, ignored = plt.hist(s, 1000, density=True)
plt.plot(bins, 1/(sigmaInc * np.sqrt(2 * np.pi)) *
           np.exp( - (bins - muInc)**2 / (2 * sigmaInc**2) ),linewidth=3, color='r')
p = np.percentile(s, 2.5) # return 50th percentile, e.g median.p = np.percentile(s, 2.5) # return 50th percentile, e.g median.
q = np.percentile(s, 97.5) # return 50th percentile, e.g median.
print(p)
print(q)

import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt

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
result = optimize.minimize(llnorm, [2,1], args = (a))
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
