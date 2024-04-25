# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 04:49:23 2024

@author: ISA
"""

import numpy as np
from setup_model import s,ref,prm,gps_born,likelihood
from obj import get_objective
from pyDOE2 import lhs
from mcmc2 import MCMC_adaptive
from tqdm import tqdm
import numpy as np
from scipy.optimize import minimize
from scipy.stats.qmc import LatinHypercube
from scipy.stats import qmc
np.random.seed(42)

def obj(x):
    return get_objective(x, ref, prm, gps_born, likelihood)

def nobj(x):
    return -obj(x)

# Flattening the bounds for Latin Hypercube Sampling
flat_bounds = []
for key, bounds in prm['bounds'].items():
    if isinstance(bounds[0], list):  # Handling nested lists for bounds
        for bound in bounds:
            flat_bounds.append(bound)
    else:
        flat_bounds.append(bounds)

# Separating lower and upper bounds
lower_bounds, upper_bounds = zip(*flat_bounds)

# Generating samples with Latin Hypercube Sampling
nsam = 1000
sampler = LatinHypercube(d=len(flat_bounds))
sample = sampler.random(n=nsam)
scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)

# Evaluating the objective function for each sample
outs = np.array([obj(x) for x in scaled_sample])

# Sorting samples by the objective function values to find the best starting point
ord_indices = np.argsort(-outs)  # Negate for descending order
best_sample = scaled_sample[ord_indices][0]

# Optimization starting from the best sample
res = minimize(nobj, best_sample, method='Nelder-Mead', options={'disp': True})

# Print the optimized parameters and the objective function value
print("Optimized Parameters:", res.x)
print("Objective Function Value:", -res.fun)