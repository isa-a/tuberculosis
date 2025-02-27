from get_objective2 import get_objective2
from setup_model import *
import numpy as np
from scipy.stats import qmc
from MCMC_adaptive import MCMC_adaptive
import matplotlib.pyplot as plt
from scipy.optimize import minimize




obj = lambda x: get_objective2(x, ref, prm, gps, prm['contmat'], lhd_fn)
nobj = lambda x: -obj(x)
nsam = 15


d = prm['bounds'].shape[1]
sampler = qmc.LatinHypercube(d=d)
lhs = sampler.random(nsam)
xsam = np.tile(prm['bounds'][0, :], (nsam, 1)) + (prm['bounds'][1, :] - prm['bounds'][0, :]) * lhs

mk = int(round(nsam / 25))
outs = np.empty(nsam)
msg = np.empty(nsam)

for ii in range(nsam):
    if (ii + 1) % mk == 0:
        print(f"{(ii + 1) / mk:.5g} ", end='')
    res, _, m_val = obj(xsam[ii, :])
    outs[ii] = res
    msg[ii] = m_val

# Order by fit (sort in descending order of outs)
mat = np.column_stack((outs, np.arange(nsam)))
sorted_idx = np.argsort(-mat[:, 0])
ord_idx = mat[sorted_idx, 1].astype(int)
xord = xsam[ord_idx, :]



# Callback function to mimic MATLAB's optimplotfval
history = []
def callback(xk):
    fxk = nobj(xk)  # Evaluate the objective (nobj) at the current iterate
    history.append(fxk)
    plt.clf()
    plt.plot(history, 'o-', markersize=4)
    plt.xlabel('Iteration')
    plt.ylabel('Objective Function Value')
    plt.title('Optimization Progress')
    plt.pause(0.01)

# Run the Nelder–Mead optimization (MATLAB's fminsearch)
# MATLAB: x0 = fminsearch(nobj, xord(1,:), options);
res = minimize(nobj, xord[0, :], method='Nelder-Mead', callback=callback, options={'disp': True})
x0 = res.x




cov0 = None
nreps = 4
niter = [int(1 * 2000), int(1 * 2000), int(1 * 2000), int(5 * 2000)]

for ii in range(nreps):
    # Call your MCMC function; note: our function returns four outputs.
    xsto, outsto, history, accept_rate = MCMC_adaptive(obj, x1, niter[ii], 1, cov0, True)
    
    # Select the parameter set corresponding to the maximum outsto
    idx = np.argmax(outsto)  # returns first index of maximum value
    x1 = xsto[idx, :]
    
    # Update the covariance matrix based on the samples (using rowvar=False so that variables are columns)
    cov0 = np.cov(xsto, rowvar=False)
    
    print()
