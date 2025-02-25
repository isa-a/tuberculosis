import numpy as np
from scipy.optimize import minimize
from scipy.stats import lognorm, beta
from scipy.special import betaln
import matplotlib.pyplot as plt

def get_distribution_fns(prctiles, distribution, visualising=False):
    dat = np.sort(np.array(prctiles))
    mn = dat[1]
    var = ((dat[2] - dat[0])**2) / 4
    target = np.array([0.025, 0.50, 0.975])
    if distribution == 'lognorm':
        def ff(x):
            mu, sigma = x
            return np.array([
                lognorm.cdf(dat[0], s=sigma, scale=np.exp(mu)),
                lognorm.cdf(dat[1], s=sigma, scale=np.exp(mu)),
                lognorm.cdf(dat[2], s=sigma, scale=np.exp(mu))
            ])
        sigma_init = np.sqrt(np.log(var / mn**2 + 1))
        mu_init = np.log(mn) - sigma_init**2 / 2
        init = [mu_init, sigma_init]
        def obj(x):
            return np.sum((ff(x) / target - 1)**2)
        res = minimize(obj, init, method='Nelder-Mead', tol=1e-12, options={'maxiter': 10**6, 'maxfev': 10**6})
        out = res.x
        val = res.fun
        if val > 1e-2:
            raise ValueError("Calibration setup not converged")
        mu_opt, sigma_opt = out
        def logfn(x):
            return -((np.log(x) - mu_opt)**2) / (2 * sigma_opt**2) - np.log(x * sigma_opt * np.sqrt(2 * np.pi))
    elif distribution == 'beta':
        def ff(x):
            a, b = x
            return np.array([
                beta.cdf(dat[0], a, b),
                beta.cdf(dat[1], a, b),
                beta.cdf(dat[2], a, b)
            ])
        tmp = mn * (1 - mn) / var - 1
        a_init = tmp * mn
        b_init = tmp * (1 - mn)
        init = [a_init, b_init]
        def obj(x):
            return np.sum((ff(x) / target - 1)**2)
        res = minimize(obj, init, method='Nelder-Mead', tol=1e-12, options={'maxiter': 10**6, 'maxfev': 10**6})
        out = res.x
        val = res.fun
        if val > 1e-2:
            raise ValueError("Calibration setup not converged")
        a_opt, b_opt = out
        def logfn(x):
            return (a_opt - 1) * np.log(x) + (b_opt - 1) * np.log(1 - x) - betaln(a_opt, b_opt)
    aux = {'sim': ff(out), 'val': val}
    if visualising:
        x_vals = np.linspace(dat[0] * 0.8, dat[2] * 1.2, 500)
        if distribution == 'lognorm':
            y_vals = lognorm.pdf(x_vals, s=out[1], scale=np.exp(out[0]))
        elif distribution == 'beta':
            y_vals = beta.pdf(x_vals, out[0], out[1])
        plt.figure()
        plt.plot(x_vals, y_vals)
        for d in dat:
            plt.axvline(x=d, linestyle='--')
        plt.show()
    return logfn, out, aux
