import numpy as np
import scipy.optimize
import logging
import fitutils as fit
__author__ = 'alex'

def conf_chi2(fun, x, y, p, nstd, sigma=1.0):
    pcopy = p
    # TODO: for some reason, this SCREWS UP the values of P given to it. Find out why.
    maxp = np.max(np.abs(p))
    conf_ints, psigma = [], []
    dbstr = ''
    def deltachi2_1dof(p0, k0, pk0):
        #chisq = fit.chi2(fun, x, y, p0, sigma=sigma)
        chisq = np.sum(np.power(y - fun(x, *p0),2))
        #pp = p0
        pp = []
        for w, pw in enumerate(p0):
            pp.append(pw)
        pp[k0] = pk0
        newfit = np.sum(np.power(y - fun(x, *pp),2))
        c2 = fit.chi2(fun, x, y, pp, sigma=sigma)
        return c2 - chisq

    for k, pk in enumerate(p):
        res = lambda dp: abs(deltachi2_1dof(pcopy, k, pk+dp) - nstd ** 2)
        #try:
        dbstr += 'before bisect: ' + str(p) + '\n'
        chi0 = scipy.optimize.minimize(res, 0, bounds=[(0,10*abs(pk))]).x[0]
        #print(chi0)
        dbstr += 'after bisect: ' + str(p) + '\n'
        #except:
        #    logging.warning('Was unable to compute a confidence limit')
        #    chi0 = 0
        conf_ints.append(2 * chi0)
        psigma.append(chi0)
    return  list(zip(p, conf_ints, psigma))

def conf1d(p, pcov, stdevs=1):
    psigma = np.sqrt(np.diagonal(pcov))
    confrad = 2*stdevs * psigma
    #print(pcov[2,:])
    r = list(zip(p, confrad, psigma))
    return r
