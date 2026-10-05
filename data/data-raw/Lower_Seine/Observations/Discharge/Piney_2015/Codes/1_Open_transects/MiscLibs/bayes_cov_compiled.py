"""bayes_cov_compiled
Computes the coefficient of variation using a Bayesian approach and an assumed posterior
log-normal distribution..

Example
-------

from MiscLibs.bayes_cov_compiled import bayes_cov

cov_68 = bayes_cov(transects, cov_prior, cov_prior_u, nsim)
"""

import numpy as np
from numba.pycc import CC
from numba import njit

cc = CC('bayes_cov_compiled')


# Bayesian COV
# ============
@cc.export('bayes_cov', 'f8(f8[::1], f8, f8, i4)')
def bayes_cov(transects_total_q, cov_prior=0.03, cov_prior_u=0.6, nsim=20000):
    """Computes the coefficient of variation using a Bayesian approach and an assumed posterior
    log-normal distribution.

    Parameters
    ----------
    transects_total_q: np.array(float)
        List of total discharge for each transect
    cov_prior: float
        Expected COV (68%) based on prior knowledge. Assumed to be 3% by default.
    cov_prior_u: float
        Uncertainty (68%) of cov_prior. Assumed to be 20%.
    nsim: int
        Number of simulations. 20000 was found to produce stable results.

    Returns
    -------
    cov: float
        Coefficient of variation
    """

    theta_std = np.abs(np.array([np.mean(transects_total_q), cov_prior])) * cov_prior_u \
        / np.sqrt(len(transects_total_q))

    # Modified for compatibility with Numba
    sam, obj_funk = metropolis(theta0=np.array([np.mean(transects_total_q), cov_prior]),
                               obs_data=transects_total_q,
                               cov_prior=cov_prior,
                               cov_prior_u=cov_prior_u,
                               nsim=nsim,
                               theta_std=theta_std)

    n_burn = int(nsim / 2)

    cov = np.mean(sam[n_burn:nsim, 1])

    return cov


@njit
@cc.export('metropolis', '(f8[:], f8[:], f8, f8, i4, f8[:])')
def metropolis(theta0, obs_data, cov_prior, cov_prior_u, nsim, theta_std):
    """Implements the Metropolis_Hastings Markov chain Monte Carlo (MCMC) algorithm for sampling the
    posterior distribution, assuming a log-normal posterior distribution.

    Parameters
    ----------
    theta0: np.array(float)
        Starting value of parameters (mean and cov_prior)
    obs_data: np.array(float)
        1D array of total discharge for each transect
    cov_prior: float
        Expected COV (68%) based on prior knowledge.
    cov_prior_u: float
        Uncertainty (68%) of cov_prior.
    nsim: int
        Number of simulations.
    theta_std: np.array(float)
        Standard deviation for the gaussian Jump distribution. If blank a default value is computed.

    Returns
    -------
    sam: np.array(float)
        Matrix containing the MCMC samples
    obj_funk: np.array(float)
        Vector containing the corresponding values of the objective function
        (i.e. of the unnormalized log-posterior)
    """

    # Initialize
    npar = len(theta0)
    sam = np.zeros((nsim + 1, npar))
    obj_funk = np.zeros((nsim + 1, 1))

    # Parameters - used for automatic computation of starting stds of the Gaussian Jump distribution
    if np.any(np.isnan(theta_std)):
        std_factor = 0.1
        theta_std = std_factor * np.abs(theta0)

    # Check if starting point is feasible - abandon otherwise
    f_current = log_post(param=theta0, measures=obs_data, cov_prior=cov_prior, cov_prior_u=cov_prior_u)

    if not is_feasible(f_current):
        print('Metropolis:FATAL:unfeasible starting point')
        return sam, obj_funk
    else:
        sam[0, :] = list(theta0)
        obj_funk[0] = f_current

        # MCMC loop
        np.random.seed(0)
        candid = np.array([np.nan, np.nan])
        for i in range(nsim):
            current = sam[i, :]
            f_current = obj_funk[i]

            # Propose a new candidate
            # Under Numba np.random.normal will not accept arrays as input
            candid[0] = np.random.normal(loc=current[0], scale=theta_std[0])
            candid[1] = np.random.normal(loc=current[1], scale=theta_std[1])

            # Evaluate objective function at candidate
            f_candid = log_post(param=candid,
                                measures=obs_data,
                                cov_prior=cov_prior,
                                cov_prior_u=cov_prior_u)

            if not is_feasible(f_candid):
                sam[i + 1, :] = current
                obj_funk[i + 1] = f_current
            else:
                # Generate deviate ~U[0,1]
                u = np.random.uniform(0, 1)

                # Compute Metropolis acceptance ratio
                ratio = np.exp(min((np.max(np.hstack((np.array([float(-100)]), f_candid - f_current))), float(0))))

                # Apply acceptance rule
                if u <= ratio:
                    sam[i + 1, :] = candid
                    obj_funk[i + 1] = f_candid
                else:
                    sam[i + 1, :] = current
                    obj_funk[i + 1] = f_current

        # Modified return for Numba, eliminating need for dictionary
        return sam, obj_funk


@njit
@cc.export('log_post', 'f8(f8[:], f8[:], f8, f8)')
def log_post(param, measures, cov_prior, cov_prior_u):
    """Define function returning the posterior log-pdf using the model measures ~ N(true_value,cov*true_value),
    with a flat prior on true_value and a log-normal prior for cov (= coefficient of variation)

    Parameters
    ----------
    param: np.array(float)
        Array containing the true value and COV

    measures: np.array(float)
        Array of observations
    cov_prior: float
        Expected COV (68%) based on prior knowledge.
    cov_prior_u: float
        Uncertainty (68%) of cov_prior.

    Returns
    -------
    logp: float
        Unnormalized log-posterior
    """
    # Check if any parameter is <=0
    # since  both true_value and cov have to be positive - otherwise sigma = true_value*cov does not make sense
    # if any(item <= 0 for item in param):
    #     return -math.inf
    # Changed for compatibility with Numba
    if np.any(np.less_equal(param, 0)):
        return np.NINF

    true_value = param[0]
    cov = param[1]
    sigma = cov * true_value  # standard deviation

    # Compute log-likelihood under the model: measures ~ N(true_value,sigma)
    # You can easily change this model (e.g. lognormal for a positive measurand?)
    # OPTION 1 : the model follows a Normal distribution
    # This equation is used for compatibility with Numba, instead of call to scipy.stats.norm.logpdf
    log_likelihood = np.sum(np.log(np.exp(-(((measures - true_value) / sigma) ** 2) / 2)
                                   / (np.sqrt(2 * np.pi) * sigma)))

    # Prior on true_value - flat prior used here but you may change this if you have prior knowledge
    log_prior_1 = 0

    # Lognormal prior
    x = cov
    mu = np.log(cov_prior)
    scale = cov_prior_u
    pdf = np.exp(-(np.log(x) - mu) ** 2 / (2 * scale ** 2)) / (x * scale * np.sqrt(2 * np.pi))
    log_prior_2 = np.log(pdf)

    # Joint prior (prior independence)
    log_prior = log_prior_1 + log_prior_2

    # Return (unnormalized) log-posterior
    logp = log_likelihood + log_prior
    if np.isnan(logp):
        # Used np to eliminate the need for math package
        logp = np.NINF  # returns -Inf rather than NaN's (required by the MCMC sampler used subsequently)
    return logp


@njit
@cc.export('is_feasible', 'b1(f8)')
def is_feasible(value):
    """Checks that a value is a real value (not infinity or nan)

    Parameters
    ----------
    value: float

    Returns
    -------
    bool
    """
    if np.isinf(value) or np.isnan(value):
        return False
    else:
        return True


if __name__ == '__main__':
    # Used to compile code
    cc.compile()
