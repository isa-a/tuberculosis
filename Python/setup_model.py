from get_addresses import get_addresses
import numpy as np
from scipy import sparse
from get_distribution_fns import get_distribution_fns
import pickle

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

# Include the auxiliaries
names = ['inc', 'incsources', 'mort', 'nTPT', 'ch_notifs']
lgths = [5, 14, 1, 1, 1]
for name, lg in zip(names, lgths):
    inds = list(range(lim + 1, lim + lg + 1))
    if 'aux' not in i:
        i['aux'] = {}
    i['aux'][name] = inds
    lim = inds[-1]
i['nx'] = lim


# -------------Set up selectors and aggregators--------------------------------

# ~~~~~~~~~Incidence
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
sel['inc'] = np.array(sel['inc'], dtype=float)



# -----------Sources of incidence 
set1 = [s['dom']]
tmp = np.zeros((len(set1), i['nstates']))
row = 0
for s_set in set1:
    inds = np.intersect1d(np.array(s['allI']), np.array(s_set))
    tmp[row, (inds - 1)] = 1
    row += 1
agg['incsources'] = sparse.csr_matrix(tmp)


# ----People starting TPT
tmp = np.zeros((i['nstates'], i['nstates']))
row_idx = np.concatenate((np.array(s['Pf']), np.array(s['Ps']))) - 1
col_idx = np.concatenate((np.array(s['Lf']), np.array(s['Ls']))) - 1
tmp[np.ix_(row_idx, col_idx)] = 1
sel['nTPT'] = sparse.csr_matrix(tmp - np.diag(np.diag(tmp)))




# -----Notifications
tmp = np.zeros((i['nstates'], i['nstates']))
sel_indices = np.intersect1d(
    np.concatenate((np.array(s['Tx']), np.array(s['Tx2']))),
    np.array(s['ch'])
) - 1
tmp[sel_indices, :] = 1
tmp[np.ix_(np.array(s['ch']) - 1, np.array(s['ad']) - 1)] = 0
tmp[np.ix_(np.array(s['ad']) - 1, np.array(s['ch']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['neg']) - 1)] = 0
tmp[np.ix_(np.array(s['neg']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['art']) - 1, np.array(s['pos']) - 1)] = 0
tmp[np.ix_(np.array(s['pos']) - 1, np.array(s['art']) - 1)] = 0
sel['ch_notifs'] = sparse.csr_matrix(tmp - np.diag(np.diag(tmp)))


# ------ Natural history parameters-------------------------------------------
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


# ---Name free parameters ----------------------------------------------------
names = ['beta', 'betadec', 'gamma', 'p_relrate_gamma_chvad', 'p_relrate', 'r_ageing_sc', 'p_relrate_factor', 'contmat_factor', 'r_ART_init', 'r_HIV_mort', 'p_HIV_relrate']
lgths = [1,         1,          2,                      1,          2,              1,          1,              1   ,                  1,           1,          1]

lim = 0
xi = {}
for name, lg in zip(names, lgths):
    inds = list(range(lim + 1, lim + lg + 1))
    xi[name] = inds
    lim = inds[-1]

bds = np.zeros((lim, 2))
bds[np.array(xi['beta']) - 1, :] = [0, 40]
bds[np.array(xi['betadec']) - 1, :] = [0, 0.15]
for idx in xi['gamma']:
    bds[idx - 1, :] = [1e-4, 10]
for idx in xi['p_relrate_gamma_chvad']:
    bds[idx - 1, :] = [0, 1]
for idx in xi['p_relrate']:
    bds[idx - 1, :] = [1, 20]
bds[np.array(xi['r_ageing_sc']) - 1, :] = [0, 1]
bds[np.array(xi['p_relrate_factor']) - 1, :] = [1, 10]
bds[np.array(xi['contmat_factor']) - 1, :] = [0, 1]
bds[np.array(xi['r_ART_init']) - 1, :] = [0, 10]
bds[np.array(xi['r_HIV_mort']) - 1, :] = [0, 10]
bds[np.array(xi['p_HIV_relrate']) - 1, :] = [1, 100]
prm = {}
prm['bounds'] = bds.T
prm['p'] = p
prm['r'] = r
prm['agg'] = agg
prm['sel'] = sel
ref = {'i': i, 's': s, 'xi': xi}

prm['contmat_born'] = np.array([[1, 0.5, 0.2],
                                [0.5, 1, 0.2],
                                [0.2, 0.2, 1]])
prm['contmat_age'] = np.array([[0.2830, 0.2525],
                               [0.0692, 0.3953]])
prm['contmat'] = np.zeros((4, 4))
for age_row in range(1, 3):
    for age_col in range(1, 3):
        for born_row in range(1, 3):
            for born_col in range(1, 3):
                row = (born_row - 1) * 2 + age_row
                col = (born_col - 1) * 2 + age_col
                prm['contmat'][row - 1, col - 1] = prm['contmat_born'][born_row - 1, born_col - 1] * prm['contmat_age'][age_row - 1, age_col - 1]


#-----------------Specify ---------------------------------------------------
data = {}
data['incd2010'] = np.array([10, 12, 14])
data['incd2020'] = np.array([6.8, 7.9, 9.2])
data['mort'] = np.array([0.26, 0.36, 0.47])
data['nTPT2019'] = 1.3 * np.array([0.9, 1, 1.1])
data['propincd_ch'] = np.array([0.006, 0.014, 0.025])
data['p_chpopn'] = np.array([0.198, 0.2471, 0.3])
data['p_adpopn'] = np.array([0.65, 0.7529, 0.85])
data['ch_notifs'] = np.array([3, 3.69, 4.2]) / 4.576e6 * 1e5
data['HIV_prev'] = np.array([0.0215, 0.0247, 0.0280])
data['ART_covg'] = np.array([0.5252, 0.6347, 0.7699])
prm['ART_start'] = 2007

ys2 = np.array([
    0, 0.0005, 0.0009, 0.0013, 0.0016, 0.0019, 0.0022, 0.0024, 0.0027, 0.0029, 0.0031, 0.0034, 0.0036, 0.0037, 0.0037, 0.0037, 0.0036, 0.0035, 0.0034,
    0.0033, 0.0031, 0.0030, 0.0028, 0.0027, 0.0026, 0.0024, 0.0023, 0.0022, 0.0021, 0.0020, 0.0019, 0.0017, 0.0016, 0.0014, 0.0013, 0.0012, 0.0012, 0.0011,
    0.0010, 0.0009, 0.0008, 0.0008, 0.0007, 0.0006, 0.0005, 0.0005, 0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
])
prm['rHIV'] = ys2

show = False
f1a = get_distribution_fns(data['incd2010'], 'lognorm', show)
f1b = get_distribution_fns(data['incd2020'], 'lognorm', show)
f2  = get_distribution_fns(data['mort'], 'lognorm', show)
f5  = get_distribution_fns(data['propincd_ch'], 'beta', show)
f6  = get_distribution_fns(data['p_chpopn'], 'beta', show)
f7  = get_distribution_fns(data['p_adpopn'], 'beta', show)
f8  = get_distribution_fns(data['ch_notifs'], 'lognorm', show)
f9  = get_distribution_fns(data['ART_covg'], 'beta', show)
f10 = get_distribution_fns(data['HIV_prev'], 'lognorm', show)


def lhd_fn(incd2010, incd2020, mort, p_chpopn, ch_notifs, ART_covg, HIV_prev):
    return (f1a(incd2010) + f1b(incd2020) + f2(mort) +
            f6(p_chpopn) + f8(ch_notifs) + f9(ART_covg) + f10(HIV_prev))


# Model_setup = {"lhd_fn": lhd_fn, "data": data, "prm": prm, "ref": ref}
# with open("Model_setup.pkl", "wb") as f:
#     pickle.dump(Model_setup, f)
