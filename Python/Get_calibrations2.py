from get_objective2 import get_objective2
from setup_model import *
import numpy as np
from scipy.stats import qmc
from make_model import make_model



obj = lambda x: get_objective2(x, ref, prm, gps, prm['contmat'], lhd)
nobj = lambda x: -obj(x)
nsam = 1000


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
