def allocate_parameters(x, p, r, prm, xi):
    import numpy as np

    r['beta'] = x[xi['beta']]
    p['betadec'] = x[xi['betadec']]
    r['gamma_2015'] = x[xi['gamma'][0]] * np.array([x[xi['p_relrate_gamma_chvad']], 1])
    r['gamma_2020'] = x[xi['gamma'][1]] * np.array([x[xi['p_relrate_gamma_chvad']], 1])

    r['progression'] = r['progression0'] * np.ones((2, 2, 3))
    mult = np.tile(np.array([[1, 0.4]]).reshape((1, 1, 2)), (2, 2, 1))
    r['progression'][:, :, 1:3] = r['progression0'] * x[xi['p_HIV_relrate']] * mult

    r['reactivation'] = r['reactivation0'] * np.ones((2, 2, 3))
    r['reactivation'][:, :, 1:3] = r['reactivation0'] * x[xi['p_HIV_relrate']] * mult

    r['ageing'] = x[xi['r_ageing_sc']]

    tmp_c = prm['contmat_age']
    prm['contmat_age'] = np.vstack((tmp_c[0, :] * x[xi['contmat_factor']], tmp_c[1, :]))

    prm['contmat'] = np.zeros((4, 4))
    for age_row in range(2):
        for age_col in range(2):
            for born_row in range(3):
                for born_col in range(3):
                    row = born_row * 2 + age_row
                    col = born_col * 2 + age_col
                    prm['contmat'][row, col] = prm['contmat_born'][born_row, born_col] * prm['contmat_age'][age_row, age_col]

    r['ART_init'] = x[xi['r_ART_init']]
    r['HIV_mort'] = x[xi['r_HIV_mort']]

    return p, r, prm
