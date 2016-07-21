__author__ = 'alex'
from math import isnan
import logging
import fitutils as fit
import scipy.optimize
import numpy as np


def shotgun_lsq(x, y, fun, p0, spread, sigma=1, maxiter = 1000):
    bestp = p0
    spread = np.array(spread)
    bestrc2 = sos(y - fun(x, *p0)) / sigma
    for k in range(maxiter):
        newp0 = ndnormal(p0, spread)
        rc2 = sos(y-fun(x,*newp0)) / sigma
        if rc2 < bestrc2:
            bestrc2 = rc2
            bestp = newp0
        if k % 500 == 0:
            spread = np.multiply(spread, 1.0)
    #print('shotgunned a reduced chi^2 of {:f}.'.format(bestrc2))
    return bestp

def ndnormal(mus, sigmas):
    """
    Creates a list of normally distributed random numbers, with each item having
    its own mean and standard deviation.
    :param n: Int
    :param mus: list of Floats
    :param sigmas: list of Floats
    :return:
    """
    r = []
    for mu, sigma in zip(mus, sigmas):
        r.append(np.random.normal(loc=mu, scale=sigma))
    return r

def sos(x):
    return np.sum(np.multiply(x, x))
