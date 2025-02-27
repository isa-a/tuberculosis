
import numpy as np
import matplotlib.pyplot as plt

def MCMC_adaptive(F, x0, n, sigma, cov0, displ):
    """
    Adaptive MCMC using Haario et al:
      F:         Function giving log-posterior density for a parameter set x
      x0:        Initial value of parameter set (1D array)
      n:         Number of iterations
      sigma:     Scaling parameter
      cov0:      Initial covariance matrix (if None, defaults to 1e-2 * I)
      displ:     Boolean flag. If True, displays progress.
      
    Returns:
      xsto:       Array of sampled parameter sets (n x d)
      outsto:     Array of log-posterior values for each sample (length n)
      history:    Array with proposed values and accept/reject flag (n x (d+1))
      accept_rate: Acceptance rate of the algorithm.
    """
    d = len(x0)
    b = 0.05
    # Compute scaling factor (MATLAB computes then overrides to 1)
    sd = sigma * (2.4**2) / d
    sd = 1  # override as in MATLAB
    
    # Check initial covariance matrix
    if cov0 is None:
        cov0 = 1e-2 * np.eye(d)
    
    # Initiate output matrices
    xsto = np.zeros((d, n))         # Each column is a sample
    outsto = np.zeros(n)
    history = np.zeros((d+1, n))      # First d rows: proposed values; last row: accept flag
    
    # Set initial state
    xsto[:, 0] = np.array(x0).flatten()
    xbar = np.copy(xsto)
    FX = F(x0)
    outsto[0] = FX
    acc = 0
    covmat_old = None
    
    if displ:
        plt.figure()
    
    # Start MCMC loop (MATLAB: for t = 2:n, here t=1,...,n-1)
    for t in range(1, n):
        X = xsto[:, t-1]
        
        # Make a proposal from the distribution:
        # Y0 ~ N(X, (0.1^2 * cov0 * sigma/d))
        try:
            Y0 = np.random.multivariate_normal(X, (0.1**2) * cov0 * sigma / d)
        except np.linalg.LinAlgError:
            Y0 = np.random.multivariate_normal(X, 1e-6 * np.eye(d))
        
        if t < 100:
            Y = np.maximum(Y0, 0)
        else:
            ind0 = t - 100
            ind1 = t - 1
            # Compute covariance from the last 100 samples:
            # Note: xsto[:, ind0:ind1+1] has shape (d, number of samples)
            covmat = np.cov(xsto[:, ind0:ind1+1])
            # Symmetrize
            covmat = (covmat + covmat.T) / 2
            
            try:
                Y = np.maximum((1 - b) * np.random.multivariate_normal(X, sd * covmat) + b * Y0, 0)
            except np.linalg.LinAlgError:
                if covmat_old is not None:
                    covmat = covmat_old
                else:
                    covmat = cov0
                Y = np.maximum((1 - b) * np.random.multivariate_normal(X, sd * covmat) + b * Y0, 0)
            covmat_old = covmat
        
        # Store proposed values in history (rows 0:d)
        history[:d, t] = Y
        
        # Decide whether to accept or reject the proposal
        FY = F(Y)
        if (np.random.rand() < np.exp(FY - FX)) and (np.abs(FY) < np.inf):
            xsel = Y
            FX = FY
            acc += 1
            history[d, t] = 1  # mark acceptance
        else:
            xsel = X
        
        xsto[:, t] = xsel
        outsto[t] = FX
        # Update running average xbar (MATLAB: xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t)
        xbar[:, t] = (xbar[:, t-1] * t + xsel) / (t + 1)
        
        # Display options
        if displ and (t % round(n / 25) == 0):
            print(f"{t/n*25:.5g} ", end='')
        if displ and (t % 200 == 0):
            plt.clf()
            for j in range(d):
                plt.plot(xsto[j, :t])
            plt.xlim([0, n])
            plt.draw()
            plt.pause(0.01)
    
    accept_rate = acc / n
    # Transpose xsto and history so that rows correspond to iterations
    xsto = xsto.T
    history = history.T
    
    return xsto, outsto, history, accept_rate
