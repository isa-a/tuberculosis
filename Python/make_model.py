import numpy as np
from scipy.sparse import csr_matrix
from setup_model import p, r, i, s, gps, prm


nstates = i['nstates']
m = np.zeros((nstates, nstates))
m2 = np.zeros((nstates, nstates))

def make_model(p, r, i, s, gps, contmat):
    
    for ia, age in enumerate(gps['age']):
        for ib, born in enumerate(gps['born']):
            for ih, hiv in enumerate(gps['hiv']):
                
                def geti(st):
                    return i[st][age][born][hiv]
                
                Lf  = geti('Lf')
                Ls  = geti('Ls')
                Pf  = geti('Pf')
                Ps  = geti('Ps')
                I   = geti('I')
                I2  = geti('I2')
                Tx  = geti('Tx')
                Tx2 = geti('Tx2')
                Rlo = geti('Rlo')
                Rhi = geti('Rhi')
                R   = geti('R')
                
                source = Lf
                destin = I
                rate = r['progression'][ia, ib, ih]
                m[destin-1, source-1] += rate
                
                source = Pf
                destin = I2
                rate = r['progression'][ia, ib, ih] * (1 - p['TPTeff'])
                m[destin-1, source-1] += rate
                
                source = Lf
                destin = Ls
                rate = r['LTBI_stabil']
                m[destin-1, source-1] += rate
                
                source = Pf
                destin = Ps
                rate = r['LTBI_stabil']
                m[destin-1, source-1] += rate
                
                source = Ls
                destin = I
                rate = r['reactivation'][ia, ib, ih]
                m[destin-1, source-1] += rate
                
                source = Ps
                destin = I
                rate = r['reactivation'][ia, ib, ih] * (1 - p['TPTeff'])
                m[destin-1, source-1] += rate
                
                source = I
                destins = [Tx, Tx2, Rhi]
                rates = [r['gamma'][ia], r['gamma'][ia], r['self_cure']]
                for d, rt in zip(destins, rates):
                    m[d-1, source-1] += rt
                
                source = I2
                destins = [Tx, Tx2, Rhi]
                rates = [r['gamma'][ia], r['gamma'][ia], r['self_cure']]
                for d, rt in zip(destins, rates):
                    m[d-1, source-1] += rt
                
                source = Tx
                destins = [Rlo, Rhi]
                rates = [r['Tx'], r['ltfu']]
                for d, rt in zip(destins, rates):
                    m[d-1, source-1] += rt
                
                source = Tx2
                destins = [Rlo, Rhi]
                rates = [r['Tx2'], r['ltfu2']]
                for d, rt in zip(destins, rates):
                    m[d-1, source-1] += rt
                
                sources = [Rlo, Rhi, R]
                destin = I2
                for src, rt in zip(sources, r['relapse']):
                    m[destin-1, src-1] += rt
                
                sources = [Rlo, Rhi]
                destin = R
                for src in sources:
                    m[destin-1, src-1] += 0.5
                
                source = Lf
                destin = Pf
                rate = r['TPT'][ib]
                m[destin-1, source-1] += rate
                
                source = Ls
                destin = Ps
                rate = r['TPT'][ib]
                m[destin-1, source-1] += rate
                
                sources = [I, I2]
                destin = Tx
                rate = r['ACF'][ib]
                for src in sources:
                    m[destin-1, src-1] += rate
                
                sources = [I, I2]
                destin = Tx2
                rate = r['ACF'][ib]
                for src in sources:
                    m[destin-1, src-1] += rate
                
                source = I2
                destin = Tx
                rate = r['ACF2'][ib]
                m[destin-1, source-1] += rate
                
                source = I2
                destin = Tx2
                rate = r['ACF2'][ib]
                m[destin-1, source-1] += rate
                
    return m

# Ageing process
sources = np.array(s['ch']) - 1
destins = np.array(s['ad']) - 1
m[destins, sources] += r['ageing']

# HIV acquisition
sources = np.array(s['neg']) - 1
destins = np.array(s['pos']) - 1
rates = 1
m2[destins, sources] += rates

# ART initiation
sources = np.array(s['pos']) - 1
destins = np.array(s['art']) - 1
rates = r['ART_init']
m[destins, sources] += rates

# Combine
M = {}
M['lin'] = csr_matrix(m - np.diag(np.sum(m, axis=0)))
M['linHIV'] = csr_matrix(m2 - np.diag(np.sum(m2, axis=0)))



# --- Nonlinear component --------------------------------------------------
M = {}
M['nlin'] = {}

for age in gps['age']:
    m = np.zeros((i['nstates'], i['nstates']))
    
    for born in gps['born']:
        for hiv in gps['hiv']:
            
            temp = np.concatenate((
                np.array(s['U']),
                np.array(s['Lf']),
                np.array(s['Ls']),
                np.array(s['Rlo']),
                np.array(s['Rhi']),
                np.array(s['R'])
            ))
            temp = np.intersect1d(temp, np.array(s[age]))
            temp = np.intersect1d(temp, np.array(s[born]))
            susinds = np.intersect1d(temp, np.array(s[hiv]))
            
            row_idx = i['Lf'][age][born][hiv] - 1
            col_indices = susinds - 1
            m[row_idx, col_indices] = 1
            
            imminds = np.concatenate((
                np.array(s['Lf']),
                np.array(s['Ls']),
                np.array(s['Rlo']),
                np.array(s['Rhi']),
                np.array(s['R'])
            )) - 1
            m[:, imminds] = m[:, imminds] * (1 - p['imm'])
            
            M['nlin'][age] = csr_matrix(m - np.diag(np.sum(m, axis=0)))
            
            
            
            
            
# --- Force of infection ---------------------------------------------------
def getinds(st1, st2):
    return np.intersect1d(
        np.intersect1d(np.array(s['infectious']), np.array(s[st1])),
        np.array(s[st2])
    )


m = np.zeros((2, i['nstates']))

m[0, getinds('ch', 'dom') - 1] = contmat[0, 0]
m[0, getinds('ad', 'dom') - 1] = contmat[0, 1]

m[1, getinds('ch', 'dom') - 1] = contmat[1, 0]
m[1, getinds('ad', 'dom') - 1] = contmat[1, 1]


# ----Include infectiousness------------------------------------------------
m = m * r['beta']
M['lam'] = csr_matrix(m)

m = np.zeros((2, i['nstates']))
m[0, np.intersect1d(np.array(s['ch']), np.array(s['dom'])) - 1] = 1
m[1, np.intersect1d(np.array(s['ad']), np.array(s['dom'])) - 1] = 1

M['denvec'] = csr_matrix(m)

# --- Mortality -----------------------------------------------------------
m = np.zeros((i['nstates'], 2))
m[np.array(s['ch']) - 1, 0] = 0
m[np.array(s['ad']) - 1, 0] = 1/72
m[np.array(s['pos']) - 1, 0] += r['HIV_mort']
m[np.array(s['infectious']) - 1, 1] = r['muTB']
M['mort'] = csr_matrix(m)
