import numpy as np
import scipy.optimize
import logging
__author__ = 'alex'

def conf_chi2(fun, x, y, p, nstd):
    # TODO: for some reason, this SCREWS UP the values of P given to it. Find out why.
    maxp = np.max(np.abs(p))
    conf_ints = []
    def deltachi2_1dof(p0, k0, pk0):
        chisq = np.sum(np.power(y - fun(x, *p0),2))
        pp = p0
        pp[k0] = pk0
        c2 = np.sum(np.power(y - fun(x, *pp),2))
        return c2 - chisq
    for k, pk in enumerate(p):
        res = lambda dp: deltachi2_1dof(p, k, pk+dp) - nstd ** 2
        try:
            chi0 = scipy.optimize.bisect(res, 0, 10*maxp)
        except:
            logging.warning('Was unable to compute a confidence limit')
            chi0 = 0
        conf_ints.append(2 * chi0)
    return  conf_ints