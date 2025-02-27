import numpy as np
import copy
from scipy.integrate import solve_ivp
from make_model2 import make_model2
from goveqs_basisnonHIV import goveqs_basisnonHIV
from goveqs_scaleup2D import goveqs_scaleup2D
from allocate_parameters import allocate_parameters
from setup_model import *


def get_objective2(x, ref, prm, gps, contmat, calfn):

    # Extract variables from ref and prm
    i = ref['i']
    s = ref['s']
    xi = ref['xi']
    p = prm['p']
    r = prm['r']
    sel = prm['sel']
    agg = prm['agg']

    # Allocate parameters (assumes allocate_parameters is defined)
    p, r, prm = allocate_parameters(x, p, r, prm, xi)

    # Check bounds conditions
    tmp1 = np.vstack((prm['bounds'], x))
    tmp2 = np.diff(tmp1[[0, 2, 1], :], axis=0)
    cond1 = np.min(tmp2) < 0

    if cond1:
        out = -np.inf
        aux = np.nan
        msg = 0
    else:
        # Set up models
        
        # Equilibrium model, pre-ART, pre-HIV
        p0 = copy.deepcopy(p)
        r0 = copy.deepcopy(r)
        p0['betadec'] = 0
        r0['gamma'] = r['gamma_2015']
        p0['relbeta'] = 0
        r0['RR_acqu'] = 0
        r0['ART_init'] = 0
        M0 = make_model2(p0, r0, i, s, gps, prm['contmat'])

        # HIV starts but no ART
        p1 = copy.deepcopy(p)
        r1 = copy.deepcopy(r)
        r1['TPT'] = [0, r['TPT2020rec'], 0, 0]
        r1['gamma'] = r['gamma_2015']
        r1['ART_init'] = 0
        M1 = make_model2(p1, r1, i, s, gps, prm['contmat'])

        # >HIV, and ART scaleup
        p2 = copy.deepcopy(p)
        r2 = copy.deepcopy(r)
        r2['gamma'] = r['gamma_2020']
        M2 = make_model2(p2, r2, i, s, gps, prm['contmat'])

        # --- Now simulate them all

        # Equilibrium model integration (pre-HIV)
        init = np.zeros(i['nx'])
        seed = 1e-5
        # Adjust indices from MATLAB (1-indexed) to Python (0-indexed)
        init[i['U']['ad']['dom']['neg'] - 1] = 1 - seed
        init[i['I']['ad']['dom']['neg'] - 1] = seed

        ode_options = {'rtol': 1e-10, 'atol': 1e-10}
        geq0 = lambda t, y: goveqs_basisnonHIV(t, y, i, s, M0, agg, sel, r0, p0, prm)
        t_span0 = (0, 5000)
        t_eval0 = np.linspace(0, 5000, 5001)
        sol0 = solve_ivp(geq0, t_span0, init, t_eval=t_eval0, **ode_options)

        # HIV decline/ART scaleup integration
        init2 = sol0.y[:, -1]
        # times: first row for gamma scaleup, second row for ART scaleup timing
        times = np.array([[2015, 2020], [2010, 2020]])
        geq1 = lambda t, y: goveqs_scaleup2D(t, y, M0, M1, M2, times, i, s, p2, r2, prm, sel, agg)
        t_span1 = (2010, 2020)
        t_eval1 = np.arange(2010, 2021)
        sol1 = solve_ivp(geq1, t_span1, init2, t_eval=t_eval1, **ode_options)

        # Compute differences along the time dimension of sol1.y
        dsol = np.diff(sol1.y, axis=1)  # shape: (nstates, len(t_eval1)-1)
        sfin = sol1.y[:, -1]

        # Record key outputs (adjust indices appropriately)
        incd2010 = dsol[i['aux']['inc'][0] - 1, 0] * 1e5
        incd2020 = dsol[i['aux']['inc'][0] - 1, -1] * 1e5
        incd = dsol[np.array(i['aux']['inc']) - 1, -1] * 1e5
        mort = dsol[np.array(i['aux']['mort']) - 1, -1] * 1e5

        n_TPT2019 = dsol[np.array(i['aux']['nTPT']) - 1, -1] * 1e5
        propincd_ch = incd[1] / incd[0] if incd[0] != 0 else np.nan
        p_chpopn = np.sum(sfin[np.array(s['ch']) - 1]) / np.sum(sfin[:i['nstates']])
        p_adpopn = np.sum(sfin[np.array(s['ad']) - 1]) / np.sum(sfin[:i['nstates']])
        ch_notifs = dsol[np.array(i['aux']['ch_notifs']) - 1, -1] * 1e5
        ART_covg = np.sum(sfin[np.array(s['art']) - 1]) / np.sum(
            sfin[np.concatenate((np.array(s['pos']), np.array(s['art']))) - 1]
        )
        HIV_prev = np.sum(sfin[np.concatenate((np.array(s['pos']), np.array(s['art']))) - 1]) / np.sum(sfin[:i['nstates']])

        if incd[0] > 0.1:
            out = calfn['fn'](incd2010, incd2020, mort, p_chpopn, ch_notifs, ART_covg, HIV_prev)
            aux = {}
            aux['soln'] = sol1.y
            msg = 2
            idx = np.where(sol1.t == 2010)[0]
            if len(idx) == 0:
                idx = 0
            else:
                idx = idx[0]
            aux['incd'] = dsol[i['aux']['inc'][0] - 1, idx:] * 1e5
            aux['incd2010'] = incd2010
            aux['incd2020'] = incd2020
            aux['mort'] = mort
            aux['nTPT'] = n_TPT2019
            aux['propincd_ch'] = propincd_ch
            aux['chpopn'] = p_chpopn
            aux['adpopn'] = p_adpopn
            aux['ch_notifs'] = ch_notifs
            aux['ART_covg'] = ART_covg
            aux['HIV_prev'] = HIV_prev
            aux['sim'] = [incd2010, incd2020, mort, propincd_ch, p_chpopn,
                          p_adpopn, ch_notifs, ART_covg, HIV_prev]
        else:
            out = -np.inf
            aux = np.nan
            msg = 1

    return out, aux, msg
