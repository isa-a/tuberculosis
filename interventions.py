# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 20:51:42 2023

@author: ISA
"""

import numpy as np
import matplotlib.pyplot as plt
from obj import get_objective
from scipy.integrate import odeint
from make_model import make_model
from allocate import allocate_parameters
from setup_model import ref, prm, sel, agg, gps_born, likelihood_function
from goveqs_basis import goveqs_basis2
from get_addresses import get_addresses


from setup_model import i,s,d,lim,r,p,agg,sel,ref,xi,prm,gps_born,likelihood
gps=gps_born
# Define the objective function obj(x)
def obj(x):
    return get_objective(x, ref, prm, gps_born, likelihood)

# Define parameters
nsam = int(1e3)
xsto = np.load('calibration_res.npy')  # Load xsto from a file
Model_setup = {}  # Define Model_setup as needed

ix0 = round(xsto.shape[0] / 2)
dx = round(xsto.shape[0] / (2 * 150))
xs = xsto[ix0::dx, :]

xx = [23.9410, 0.1473, 5.2634, 0.7268, 12.3779]

mk = round(xs.shape[0] / 25)
outs = np.zeros(nsam)

out, aux = obj(xx)
init = aux['soln'][-1, :]

# Allocate parameters (replace with actual implementation)
p0, r0 = allocate_parameters(xx, p, r, xi)

# Make the model (replace with actual implementation)
M0 = make_model(p0, r0, i, s, gps)

# Define other models (replace with actual implementations)
M1 = make_model(p1, r1, i, s, gps)
M2 = make_model(p2, r2, i, s, gps)
M3 = make_model(p3, r3, i, s, gps)
M4 = make_model(p4, r4, i, s, gps)
M5 = make_model(p5, r5, i, s, gps)

# Create a list of models
models = [M0, M1, M2, M3, M4, M5]

for mi in range(len(models)):
    # Define the geq function (replace with actual implementation)
    def geq(t, in_values):
        return goveqs_scaleup(t, in_values, i, s, M0, models[mi], p0, p1, [2022, 2025], agg, sel, r)

    # Solve the ODE (replace with actual solver)
    t, soln = ode_solver(geq, [2022, 2036], init)

    # Calculate and store results (replace with actual calculations)
    sdiff = np.diff(soln, axis=0)
    incsto[:, ii, mi] = sdiff[:, i['aux']['inc'][0]] * 1e5
    mrtsto[:, ii, mi] = sdiff[:, i['aux']['mort'][0]] * 1e5

    # Calculate proportions (replace with actual calculations)
    vec = sdiff[-1, i['aux']['incsources']] * 1e5
    props[ii, :, mi] = vec / np.sum(vec)

print('\n')

# Calculate percentiles for incidence and mortality
incmat = np.percentile(incsto, [2.5, 50, 97.5], axis=1).transpose((1, 0, 2))
mrtmat = np.percentile(mrtsto, [2.5, 50, 97.5], axis=1).transpose((1, 0, 2))

# Plot incidence
cols = plt.cm.viridis(np.linspace(0, 1, incmat.shape[2]))
plt.figure()
lw = 1.5
fs = 14
xx = np.arange(2022, 2036)
for ii in range(incmat.shape[2]):
    plt.plot(xx, incmat[:, 1, ii], color=cols[ii], linewidth=lw)
    plt.fill_between(xx, incmat[:, 2, ii], incmat[:, 0, ii], color=cols[ii], alpha=0.1)
    plt.axhline(y=1.05, color='k', linestyle='--')
yl = plt.ylim()
yl = (0, yl[1])
plt.ylim(yl)
plt.xlim((2022, 2035))
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend(['Baseline', 'ACF', 'ACF + domestic TPT', 'ACF + domestic AND migrant TPT', '+ Monthly followup post TPT or Tx', 'Elimination target'], loc='southwest')
plt.ylabel('Rate per 100,000 population')
plt.title('Incidence')
plt.show()

# Plot mortality
cols = plt.cm.viridis(np.linspace(0, 1, mrtmat.shape[2]))
plt.figure()
for ii in range(mrtmat.shape[2]):
    plt.plot(xx, mrtmat[:, 1, ii], color=cols[ii], linewidth=lw)
    plt.fill_between(xx, mrtmat[:, 2, ii], mrtmat[:, 0, ii], color=cols[ii], alpha=0.1)
yl = plt.ylim()
yl = (0, yl[1])
plt.ylim(yl)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend(['Baseline', 'ACF', 'ACF + domestic TPT', 'ACF + domestic AND migrant TPT', '+ Monthly followup post TPT or Tx', 'Elimination target'], loc='southwest')
plt.ylabel('Rate per 100,000 population')
plt.title('Mortality')
plt.show()

# Show proportions from different sources
tmp1 = np.percentile(props, [2.5, 50, 97.5], axis=0)
tmp2 = tmp1[1, :, -1]
mm = [np.sum(tmp2[[0, 2]]), np.sum(tmp2[[1, 3]]), tmp2[4]]

plt.figure()
plt.pie(mm, labels=['People developing TB without history of TPT or active TB treatment', 'People developing TB after TPT', 'Relapse after active TB treatment'])
plt.legend(loc='northwest', orientation='vertical')
plt.title('Sources of incidence in 2035 with all interventions combined')
plt.show()
