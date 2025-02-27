import numpy as np
from scipy.sparse import csr_matrix
from setup_model import p, r, i, s, gps, prm


def make_model2(p, r, i, s, gps, contmat):
    import numpy as np
    from scipy.sparse import csr_matrix

    nstates = i['nstates']
    m = np.zeros((nstates, nstates))
    m2 = np.zeros((nstates, nstates))

    # --- Linear transitions: Loop over age, born, HIV status
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
                
                # Progression from 'fast' latent
                source = Lf
                destin = I
                rate = r['progression'][ia, ib, ih]
                m[destin - 1, source - 1] += rate
                
                source = Pf
                destin = I2
                rate = r['progression'][ia, ib, ih] * (1 - p['TPTeff'])
                m[destin - 1, source - 1] += rate
                
                # Stabilisation of 'fast' to 'slow' latent
                source = Lf
                destin = Ls
                rate = r['LTBI_stabil']
                m[destin - 1, source - 1] += rate
                
                source = Pf
                destin = Ps
                rate = r['LTBI_stabil']
                m[destin - 1, source - 1] += rate
                
                # Reactivation of 'slow' latent
                source = Ls
                destin = I
                rate = r['reactivation'][ia, ib, ih]
                m[destin - 1, source - 1] += rate
                
                source = Ps
                destin = I
                rate = r['reactivation'][ia, ib, ih] * (1 - p['TPTeff'])
                m[destin - 1, source - 1] += rate
                
                # Initiation of treatment
                source = I
                destins = [Tx, Tx2, Rhi]
                rates = [r['gamma'][ia], r['gamma'][ia], r['self_cure']]
                for d, rt in zip(destins, rates):
                    m[d - 1, source - 1] += rt
                
                source = I2
                destins = [Tx, Tx2, Rhi]
                rates = [r['gamma'][ia], r['gamma'][ia], r['self_cure']]
                for d, rt in zip(destins, rates):
                    m[d - 1, source - 1] += rt
                
                # Treatment completion or interruption
                source = Tx
                destins = [Rlo, Rhi]
                rates = [r['Tx'], r['ltfu']]
                for d, rt in zip(destins, rates):
                    m[d - 1, source - 1] += rt
                
                # Second-line treatment
                source = Tx2
                destins = [Rlo, Rhi]
                rates = [r['Tx2'], r['ltfu2']]
                for d, rt in zip(destins, rates):
                    m[d - 1, source - 1] += rt
                
                # Relapse
                sources_arr = [Rlo, Rhi, R]
                destin = I2
                for src, rt in zip(sources_arr, r['relapse']):
                    m[destin - 1, src - 1] += rt
                
                # Stabilisation of relapse risk
                sources_arr = [Rlo, Rhi]
                destin = R
                for src in sources_arr:
                    m[destin - 1, src - 1] += 0.5
                
                # Initiation of TPT
                source = Lf
                destin = Pf
                rate = r['TPT'][ib]
                m[destin - 1, source - 1] += rate
                
                source = Ls
                destin = Ps
                rate = r['TPT'][ib]
                m[destin - 1, source - 1] += rate
                
                # Case-finding
                sources_arr = [I, I2]
                destin = Tx
                rate = r['ACF'][ib]
                for src in sources_arr:
                    m[destin - 1, src - 1] += rate
                
                sources_arr = [I, I2]
                destin = Tx2
                rate = r['ACF'][ib]
                for src in sources_arr:
                    m[destin - 1, src - 1] += rate
                
                # ACF2
                source = I2
                destin = Tx
                rate = r['ACF2'][ib]
                m[destin - 1, source - 1] += rate
                
                source = I2
                destin = Tx2
                rate = r['ACF2'][ib]
                m[destin - 1, source - 1] += rate

    # --- Ageing process
    sources = np.array(s['ch']) - 1
    destins = np.array(s['ad']) - 1
    m[np.ix_(destins, sources)] += r.get('ageing', 0)
    
    # --- HIV acquisition
    sources = np.array(s['neg']) - 1
    destins = np.array(s['pos']) - 1
    m2[np.ix_(destins, sources)] += 1
    
    # --- ART initiation
    sources = np.array(s['pos']) - 1
    destins = np.array(s['art']) - 1
    m[np.ix_(destins, sources)] += r['ART_init']
    
    # --- Combine linear components
    M_lin = csr_matrix(m - np.diag(np.sum(m, axis=0)))
    M_linHIV = csr_matrix(m2 - np.diag(np.sum(m2, axis=0)))
    
    # --- Nonlinear component
    M_nlin = {}
    for age in gps['age']:
        m_nlin = np.zeros((nstates, nstates))
        for born in gps['born']:
            for hiv in gps['hiv']:
                arr1 = np.concatenate((
                    np.array(s['U']),
                    np.array(s['Lf']),
                    np.array(s['Ls']),
                    np.array(s['Rlo']),
                    np.array(s['Rhi']),
                    np.array(s['R'])
                ))
                arr2 = np.array(s[age])
                arr3 = np.array(s[born])
                arr4 = np.array(s[hiv])
                temp = np.intersect1d(arr1, arr2)
                temp = np.intersect1d(temp, arr3)
                susinds = np.intersect1d(temp, arr4)
                
                row_idx = i['Lf'][age][born][hiv] - 1
                col_indices = susinds - 1
                m_nlin[row_idx, col_indices] = 1
                
                imminds = np.concatenate((
                    np.array(s['Lf']),
                    np.array(s['Ls']),
                    np.array(s['Rlo']),
                    np.array(s['Rhi']),
                    np.array(s['R'])
                )) - 1
                m_nlin[:, imminds] = m_nlin[:, imminds] * (1 - p['imm'])
        M_nlin[age] = csr_matrix(m_nlin - np.diag(np.sum(m_nlin, axis=0)))
    
    # --- Force of infection
    def getinds(st1, st2):
        return np.intersect1d(
            np.intersect1d(np.array(s['infectious']), np.array(s[st1])),
            np.array(s[st2])
        )
    
    # Ensure contmat exists. The MATLAB code sets contmat(end,end)=contmat(end,end)
    contmat[-1, -1] = contmat[-1, -1]
    
    m_foi = np.zeros((2, nstates))
    m_foi[0, getinds('ch', 'dom') - 1] = contmat[0, 0]
    m_foi[0, getinds('ad', 'dom') - 1] = contmat[0, 1]
    m_foi[1, getinds('ch', 'dom') - 1] = contmat[1, 0]
    m_foi[1, getinds('ad', 'dom') - 1] = contmat[1, 1]
    
    # --- Include infectiousness
    m_foi = m_foi * r['beta']
    M_lam = csr_matrix(m_foi)
    
    # --- Additional matrix to track numbers in each group (denvec)
    m_den = np.zeros((2, nstates))
    m_den[0, np.intersect1d(np.array(s['ch']), np.array(s['dom'])) - 1] = 1
    m_den[1, np.intersect1d(np.array(s['ad']), np.array(s['dom'])) - 1] = 1
    M_denvec = csr_matrix(m_den)
    
    # --- Mortality
    m_mort = np.zeros((nstates, 2))
    m_mort[np.array(s['ch']) - 1, 0] = 0
    m_mort[np.array(s['ad']) - 1, 0] = 1/72
    m_mort[np.array(s['pos']) - 1, 0] += r['HIV_mort']
    m_mort[np.array(s['infectious']) - 1, 1] = r['muTB']
    M_mort = csr_matrix(m_mort)
    
    # --- Combine all components into output dictionary M
    M = {}
    M['lin'] = M_lin
    M['linHIV'] = M_linHIV
    M['nlin'] = M_nlin
    M['lam'] = M_lam
    M['denvec'] = M_denvec
    M['mort'] = M_mort

    return M
