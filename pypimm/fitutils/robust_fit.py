__author__ = 'alex'
import scipy
import scipy.optimize
import numpy as np
from math import isnan, pi
import logging
import matplotlib.pyplot as plt

def robust_fit(x, y, fun, p0, spread, dist=None, sigma=1, maxiter=1000):
    """
    Fits data to a model using robust initial guess and
    parameter estimation tecnhiques.

    Initial guess refinement is done with shotgun least-squares,
    and model fitting is done via maximum-likelihood estimates
    powered by Nelder-Mead cumulative error minimization.

    :param x:
    :param y:
    :param fun:
    :param p0:
    :param spread:
    :param dist:
    :param sigma:
    :param maxiter:
    :return:
    """
    if not dist:
        dist = lambda xi: np.multiply(xi, xi)
    cumres = lambda p: np.sum(dist(y - fun(x, *p)))

    bestguess = shotgun_lsq(x, y, fun, p0, spread, dist=dist, sigma=sigma, maxiter=maxiter)
    optres = scipy.optimize.minimize(cumres, bestguess, method='Nelder-Mead')
    print(optres)
    bestp = optres.x[0]
    bestfit = fun(x, *bestp)
    ssres = sum(dist(y - bestfit))  # sum of squares of the residuals

    ymean = np.mean(y)
    sstot = sum(dist(y - ymean) )   # total sum of squares
    r = 1.0 - ssres / sstot
    bestcov = 0
    return bestp, bestcov, r

def shotgun_lsq(x, y, fun, p0, spread, sigma=1, dist=None, maxiter = 1000):
    if not dist:
        cumres = lambda x: np.sum(np.multiply(x, x))
    else:
        cumres = lambda x: np.sum(dist(x))
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
            mse = cumres(y - bestfit) / sigma
        except ValueError:
            mse = 1E9
        if(isnan(mse)):
            mse = 1E9
        if mse < bestmse:
            bestmse = mse
            bestp = newp0
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

def dsinplus_sp(x, p0, p1, p2, p3, p4, p5):
    """
    Exactly the same as dsin, but with separate args because that's what Scipy's curve_fit takes.
    :param x:
    :param p:
    :return:
    """
    x = np.array(x)
    y =   p0 * np.cos((x-p5)*2*pi*p1)*np.exp(-(x-p5)/abs(p2)) \
        + p3 * np.exp(-(x-p5)/abs(p4))
    y[y == np.inf] = 0
    y[y == np.nan] = 0
    return y

def ml_lorentz(z):
    return np.log(1 + 0.5 * z * z)


t = np.linspace(0, 5, 512)
y = dsinplus_sp(t, 2, 3, 1, -3, 1, 0) + np.random.random(512)
p0 = [1,1,1,1,1,1]
spr = [1,1,1,1,1,1]
pestimate = robust_fit(t, y, dsinplus_sp, p0, spread=spr, dist=None)
print(pestimate)
yfit = dsinplus_sp(t, *pestimate)
plt.plot(t,y, t, yfit)
plt.show()