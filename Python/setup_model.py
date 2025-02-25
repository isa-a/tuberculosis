from get_addresses import get_addresses
import numpy as np
from scipy import sparse


states1 = ['U']
states2 = ['Lf', 'Ls', 'Pf', 'Ps', 'I', 'I2', 'Tx', 'Tx2', 'Rlo', 'Rhi', 'R']

gps = {
    'age': ['ch', 'ad'],
    'born': ['dom'],
    'hiv': ['neg', 'pos', 'art']
}

i, s, d, lim = get_addresses([states1, gps['age'], gps['born'], gps['hiv']], {}, {}, {}, 0)
i, s, d, lim = get_addresses([states2, gps['age'], gps['born'], gps['hiv']], i, s, d, lim)
d = {k: ' '.join(v) for k, v in d.items()}

s['allI'] = s['I'] + s['I2']
s['infectious'] = s['allI'] + s['Tx']
s['prev'] = s['allI'] + s['Tx']

names = ['inc', 'incsources', 'mort', 'nTPT', 'ch_notifs']
lgths = [5, 14, 1, 1, 1]
for name, lg in zip(names, lgths):
    inds = list(range(lim + 1, lim + lg + 1))
    if 'aux' not in i:
        i['aux'] = {}
    i['aux'][name] = inds
    lim = inds[-1]
i['nx'] = lim


tmp = np.zeros((5, i['nstates']))
tmp[0, np.array(s['allI']) - 1] = 1
tmp[1, np.array(np.intersect1d(s['allI'], s['ch'])) - 1] = 1
tmp[2, np.array(np.intersect1d(s['allI'], s['neg'])) - 1] = 1
tmp[3, np.array(np.intersect1d(s['allI'], s['pos'])) - 1] = 1
tmp[4, np.array(np.intersect1d(s['allI'], s['art'])) - 1] = 1
agg = {}
agg['inc'] = sparse.csr_matrix(tmp)

tmp = np.zeros((i['nstates'], i['nstates']))
tmp[np.array(s['allI']) - 1, :] = 1
tmp[np.ix_(np.array(s['ch']) - 1, np.array(s['ad']) - 1)] = 0
tmp[np.ix_(np.array(s['ad']) - 1, np.array(s['ch']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['neg']) - 1)] = 0
tmp[np.ix_(np.array(s['neg']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['art']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['art']) - 1)] = 0
sel = {}
sel['inc'] = tmp - np.diag(np.diag(tmp))

set1 = [s['dom']]
tmp = np.zeros((len(set1), i['nstates']))
row = 0
for s_set in set1:
    inds = np.intersect1d(np.array(s['allI']), np.array(s_set))
    tmp[row, (inds - 1)] = 1
    row += 1
agg['incsources'] = sparse.csr_matrix(tmp)

tmp = np.zeros((i['nstates'], i['nstates']))
row_idx = np.concatenate((np.array(s['Pf']), np.array(s['Ps']))) - 1
col_idx = np.concatenate((np.array(s['Lf']), np.array(s['Ls']))) - 1
tmp[np.ix_(row_idx, col_idx)] = 1
sel['nTPT'] = tmp - np.diag(np.diag(tmp))

tmp = np.zeros((i['nstates'], i['nstates']))
sel_indices = np.intersect1d(np.concatenate((np.array(s['Tx']), np.array(s['Tx2']))), np.array(s['ch'])) - 1
tmp[sel_indices, :] = 1
tmp[np.ix_(np.array(s['ch']) - 1, np.array(s['ad']) - 1)] = 0
tmp[np.ix_(np.array(s['ad']) - 1, np.array(s['ch']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['neg']) - 1)] = 0
tmp[np.ix_(np.array(s['neg']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['art']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['art']) - 1)] = 0
sel['ch_notifs'] = tmp - np.diag(np.diag(tmp))

r = {}
p = {}

r['progression0'] = 0.0826
r['LTBI_stabil'] = 0.872
r['reactivation0'] = 0.0006
r['Tx'] = 2
p['RRrec'] = 1
r['RR_acqu'] = 0
r['Tx2'] = 9/12
r['ltfu'] = 0.01
r['ltfu2'] = r['Tx2'] * 2
r['self_cure'] = 1/6
r['relapse'] = [0.032, 0.14, 0.0015]
r['muTB'] = 1/6
p['imm'] = 0.8
p['migrTPT'] = 0
p['TPTeff'] = 0.6
r['TPT'] = [0, 0, 0, 0]
r['TPT2020rec'] = 0.004
r['ACF'] = [0, 0, 0, 0]
r['ACF2'] = [0, 0, 0, 0]
