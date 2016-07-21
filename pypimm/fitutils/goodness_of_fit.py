import numpy as np
import scipy.optimize
__author__ = 'alex'

def gof(fun, x, y, p):
    """
    Computes a goodness-of-fit value for a given model estimate
    to a set of data. Currently, this consists of reduced chi squared.
    :param fun:
    :param y:
    :param p:
    :return:
    """

    return 1

def const_chi(fun, x, y, p, conf_int=0.95):
    """
    compute the constant-chi confidence interval for each
    individual parameter in a model fit for a set of data.
    :param fun:
    :param x:
    :param y:
    :param p:
    :param conf_int:
    :return:
    """

    return None

def minimize_reduced_chi2(fun, x, y, p0, sigma=1.0):
    res = lambda p: abs(np.add(redchi2(fun, x, y, p, sigma=sigma), -1.0))
    popt = scipy.optimize.minimize(res, p0, method='Nelder-Mead').x
    print(popt)
    bestchi2 = res(popt) + 1
    print(bestchi2)
    return popt, bestchi2

def minimize_absolute(fun, x, y, p0, sigma=1.0):
    res = lambda p: np.sum(np.abs(y - fun(x, *p))/sigma)
    popt = scipy.optimize.minimize(res, p0, method='Nelder-Mead').x
    bestchi2 = chi2(fun, x, y, popt, sigma=sigma)
    print(popt)
    return popt, bestchi2

def basin_lsq(fun, x, y, p0, sigma=1.0):
    res = lambda p: np.sum(np.power((y - fun(x, *p)),2))
    popt = scipy.optimize.basinhopping(res, p0, niter=150).x
    bestchi2 = chi2(fun, x, y, popt, sigma=sigma)
    return popt, bestchi2

def minimize_lorentz(fun, x, y, p0, sigma=1.0):
    def res(p):
        a = np.power((y - fun(x, *p))/sigma, 2)
        return np.sum(np.log(0.5*a + 1))
    popt = scipy.optimize.minimize(res, p0, method='Nelder-Mead').x
    bestchi2 = chi2(fun, x, y, popt, sigma=sigma)
    print(popt)
    return popt, bestchi2

def redchi2(fun, x, y, p, sigma=1.0):
    n = np.size(y)
    m = np.size (p)
    a = np.divide(y - fun(x, *p), sigma)
    return np.sum(np.power(a,2)) / (n - m)

def chi2(fun, x, y, p, sigma=1.0):
    a = np.divide(y - fun(x, *p), sigma)
    return np.sum(np.power(a,2))

def nlcorr(fun, x, y, p, sigma=1.0):
    bestfit = fun(x, *p)
    ssres = np.sum(np.power(y - bestfit, 2))  # sum of squares of the residuals
    sigmean = np.mean(y)
    sstot = np.sum(np.power(y - sigmean, 2))    # total sum of squares
    r2 = 1.0 - ssres / sstot
    return r2