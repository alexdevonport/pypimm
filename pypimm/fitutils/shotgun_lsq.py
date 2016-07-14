__author__ = 'alex'
from math import isnan
import logging

import scipy.optimize
import numpy as np


def shotgun_lsq(x, y, fun, p0, spread, sigma=1, maxiter = 1000):
    mse = 1
    bestmse = 2
    m = len(p0)
    bestp = p0
    yrange = np.max(y) - np.min(y)
    spread = np.array(spread)
    for k in range(maxiter):
        newp0 = ndnormal(p0, spread)
        try:
            bestfit = fun(x, *newp0)
            mse = sos(y - bestfit) / sigma
        except ValueError:
            mse = 1E9
        if(isnan(mse)):
            mse = 1E9
        if mse < bestmse:
            bestmse = mse
            bestp = newp0
        if k % 250 == 0:
            spread = np.multiply(spread, 1.25)
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
