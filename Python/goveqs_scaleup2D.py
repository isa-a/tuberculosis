from goveqs_basis2 import goveqs_basis2
import numpy as np

def goveqs_scaleup2D(t, in_vec, M0, M1, M2, times, i, s, p, r, prm, sel, agg):

    scale = np.maximum((t - times[:, 0]) / (times[:, 1] - times[:, 0]), 0)
    scale[0] = min(scale[0], 1)
    Mt = M1.copy()
    Mt['lin'] = M0['lin'] + scale[0]*(M1['lin'] - M0['lin']) + scale[1]*(M2['lin'] - M0['lin'])
    out = goveqs_basis2(t, in_vec, i, s, Mt, agg, sel, r, p, prm)
    return out
