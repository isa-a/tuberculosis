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
#import model
from pymcmcstat.MCMC import MCMC

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


solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
U, Lf, Ls, I, R, cInc = solve.T #get trajectories


def loglik(beta, gamma):
    
    solve = odeint(deriv, (U0, Lf0, Ls0, I0, R0, J0), t, args=(N, beta, gamma, mu, muTB, sigma, rho, u, v, w))
    U, Lf, Ls, I, R, cInc = solve.T #get trajectories

    muPrev, sigmaPrev = I[-1]*100000, 40 #I (prevalence)
    muInc, sigmaInc = (cInc[1:] - cInc[:-1])[-1]*100000, 30 #cInc (incidence)
    n = 10000
    
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



mcstat = MCMC()  #initialise mcmc
x = t
y = I*100000
mcstat.data.add_data_set(x,y) #create data using prevalence



mcstat.model_settings.define_model_settings(sos_function=loglik)  # put function into model
mcstat.simulation_options.define_simulation_options(nsimu=10.0e3) #number of sims to run
mcstat.parameters.add_model_parameter(name='beta', theta0=8, minimum = 0, maximum = 20) #priors
mcstat.parameters.add_model_parameter(name='gamma',theta0=0.4, minimum = 0, maximum=2) #priors
mcstat.run_simulation() #run mcmc

results = mcstat.simulation_results.results #get results
chain = results['chain'] #get results
burnin = int(chain.shape[0]/2) #define burnin period as first 50% of chain
# display chain statistics
mcstat.chainstats(chain[burnin:, :], results)

names = results['names']
chain = results['chain']
s2chain = results['s2chain']
names = results['names'] # parameter names

mcpl = mcstat.mcmcplot # initialize plotting methods
mcpl.plot_chain_panel(chain, names)

mcpl.plot_chain_panel(chain[burnin:,:], names)

nds = 100
x = np . linspace (2 , 3 , num = nds )
x = x . reshape ( nds ,1)
m = 2 # slope
b = -3 # offset
noise = 0.1* np . random . standard_normal ( x . shape )
y = m * x + b + noise
mcstat = MCMC ()
mcstat.data.add_data_set(x, y)
mcstat.parameters.add_model_parameter(name = 'm', theta0 = 1. , minimum = -10 , maximum = 10)
mcstat.parameters.add_model_parameter ( name = 'b', theta0 = -5. , minimum =-10 , maximum = 100)
mcstat.simulation_options.define_simulation_options(nsimu=10.0e3) #number of sims to run
def test_modelfun ( xdata , theta ) :
    m = theta [0]
    b = theta [1]
    nrow , ncol = xdata . shape
    y = np . zeros ([ nrow ,1])
    y [: ,0] = m * xdata . reshape ( nrow ,) + b
    return y
def test_ssfun ( theta , data ) :
    xdata = data . xdata [0]
    ydata = data . ydata [0]
    # eval model
    ymodel = test_modelfun ( xdata , theta )
    # calc sos
    ss = sum (( ymodel [: ,0] - ydata [: ,0]) **2)
    return ss
####################### NIM: IGNORE BELOW ###########################################


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
###################################################################################################################################################

import numpy as np
import scipy
import scipy.stats
import matplotlib as mpl   
import matplotlib.pyplot as plt

mod1=lambda t:np.random.normal(10,3,t)

#Form a population of 30,000 individual, with average=10 and scale=3
population = mod1(30000)
#Assume we are only able to observe 1,000 of these individuals.
observation = np.random.normal(400, 40, 1000)
 

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.hist( observation,bins=35 ,)
ax.set_xlabel("Value")
ax.set_ylabel("Frequency")
ax.set_title("Figure 1: Distribution of 1000 observations sampled from a population of 30,000 with $\mu$=10, $\sigma$=3")
mu_obs=observation.mean()
mu_obs

transition_model = lambda x: [x[0],np.random.normal(x[1],0.5,(1,))[0]]

def prior(x):
    #x[0] = mu, x[1]=sigma (new or current)
    #returns 1 for all valid values of sigma. Log(1) =0, so it does not affect the summation.
    #returns 0 for all invalid values of sigma (<=0). Log(0)=-infinity, and Log(negative number) is undefined.
    #It makes the new sigma infinitely unlikely.
    if(x[1] <=0):
        return 0
    return 1

#Computes the likelihood of the data given a sigma (new or current) according to equation (2)
def manual_log_like_normal(x,data):
    #x[0]=mu, x[1]=sigma (new or current)
    #data = the observation
    return np.sum(-np.log(x[1] * np.sqrt(2* np.pi) )-((data-x[0])**2) / (2*x[1]**2))


#Defines whether to accept or reject the new sample
def acceptance(x, x_new):
    if x_new>x:
        return True
    else:
        accept=np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        # less likely x_new are less likely to be accepted
        return (accept < (np.exp(x_new-x)))


def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,data,acceptance_rule):
    # likelihood_computer(x,data): returns the likelihood that these parameters generated the data
    # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(x,x_new): decides whether to accept or reject the new sample
    x = param_init
    accepted = []
    rejected = []   
    for i in range(iterations):
        x_new =  transition_model(x)    
        x_lik = likelihood_computer(x,data)
        x_new_lik = likelihood_computer(x_new,data) 
        if (acceptance_rule(x_lik + np.log(prior(x)),x_new_lik+np.log(prior(x_new)))):            
            x = x_new
            accepted.append(x_new)
        else:
            rejected.append(x_new)            
                
    return np.array(accepted), np.array(rejected)
accepted, rejected = metropolis_hastings(manual_log_like_normal,prior,transition_model,[mu_obs,0.1], 50000,observation,acceptance)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2,1,1)

ax.plot( rejected[0:100,1], 'rx', label='Rejected',alpha=0.5)
ax.plot( accepted[0:100,1], 'b.', label='Accepted',alpha=0.5)
ax.set_xlabel("Iteration")
ax.set_ylabel("$\sigma$")
ax.set_title("Figure 2: MCMC sampling for $\sigma$ with Metropolis-Hastings. First 50 samples are shown.")
ax.grid()
ax.legend()



ax2 = fig.add_subplot(2,1,2)
to_show=-accepted.shape[0]
ax2.plot( rejected[to_show:,1], 'rx', label='Rejected',alpha=0.5)
ax2.plot( accepted[to_show:,1], 'b.', label='Accepted',alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("$\sigma$")
ax2.set_title("Figure 3: MCMC sampling for $\sigma$ with Metropolis-Hastings. All samples are shown.")
ax2.grid()
ax2.legend()



fig.tight_layout()
accepted.shape


show=int(-0.75*accepted.shape[0])
hist_show=int(-0.75*accepted.shape[0])

fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(1,2,1)
ax.plot(accepted[show:,1])
ax.set_title("Figure 4: Trace for $\sigma$")
ax.set_ylabel("$\sigma$")
ax.set_xlabel("Iteration")
ax = fig.add_subplot(1,2,2)
ax.hist(accepted[hist_show:,1], bins=20,density=True)
ax.set_ylabel("Frequency (normed)")
ax.set_xlabel("$\sigma$")
ax.set_title("Figure 5: Histogram of $\sigma$")
fig.tight_layout()


ax.grid("off")
