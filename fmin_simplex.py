__author__ = 'alex'
import numpy as np
import logging
from math import pi, isnan
import matplotlib.pyplot as plt
import scipy.optimize

logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s - %(message)s')

def fmin_simplex(fun, p0, spread=1, rtol = 1E-10):
    """
    Finds a local minimum of
    :param fun: function that takes a list of length M floats as an argument
                and returns a float
    :param p0: list of length M
    :param spread: float
    :param rtol: float
    :return:tuple of (float, [float, float])
    """
    d = len(p0)
    npts = d + 1
    # initialize simplex and auxiliary variables.
    # first index is points on the simplex, so p[2,:]
    # is point 2 on the simplex.
    p = np.zeros((npts,d))
    dist = 0
    fs = np.zeros(npts)
    tryno = 0
    for row in p:
        row += p0 + spread * np.random.rand(d)
    # The matrix p now represents a simplex with points around p0.
    satisfied = False   # We will be satisfied when we find a minimum
    lastmin = 0
    stationarycount = 0
    while not satisfied:
        debugstr = ''
        tryno += 1
        debugstr += ' Iteration #'+str(tryno)+'. '
        # compute the best and worst points
        for k,row in enumerate(p):
            fs[k] = fun(row) # find function value for each point
        current_best = np.min(fs)
        best, worst = np.argmin(fs), np.argmax(fs)  # indices of the best and worst points
        current_worst = p[worst] # used to calculate distance traveled this step
        # calculate the centroid of all points but the worst
        m = np.zeros(d)
        for row in np.delete(p, worst, axis = 0): # iterates over all rows but the worst
            m += row / (npts - 1)
        # with the centroid, we can calculate the 'reflection' of the worst point
        # and its function value
        r = 2 * m - p[worst]
        fr = fun(r)
        if fs[worst] > fr > fs[best]:
            # Case 1: the new point is neither the best nor the worst.
            # Replace the worst point with the new one. It's an improvement, at least.
            dist = abs(r - p[worst])
            p[worst] = r
            debugstr += ' Attempting reflection.'
        elif fr < fs[best]:
            # Case 2: the new point is better than all of them. This is probably
            # the right direction! We'll make an expanded point (which is 'reflected'
            # twice the distance from the centroid). If that's better than the original
            # reflected point, we'll replace the worst one with that; otherwise, we'll
            # just stick with the original.
            e = 2 * r - m
            if fun(e) < fr:
                p[worst] = e
                debugstr += ' Attempting expansion.'
            else:
                p[worst] = r
                debugstr += ' Attempted expansion, but failed. Reflecting instead.'
        else:
            # It looks like the new point is at least worst than the second
            # worst point. We'll try contracting the reflected point, to see
            # if that's any better. If it isn't, the last thing to do is to shrink
            # all points.
            c = 0.5 * m + 0.5 * p[worst]
            if fun(c) < fs[worst]:
                p[worst] = c
                debugstr += ' Attempting contraction.'
            else:
                debugstr += ' Attempting shrink.'
                newp = np.zeros(np.shape(p))
                for k, row in enumerate(p):
                    if k != best:
                        newp[k] = 0.5 * p[best] + 0.5 * row
                p = newp
        # Now that we have reflected, expanded, contracted, or shrunk, we can see
        # what the new best value is, and if we can stop.
        for k,row in enumerate(p):
            fs[k] = fun(row) # find function value for each point
        fmin = np.min(fs)
        #if abs(fmin - lastmin) < rtol:
        #    print('STATIONARY')
        #    stationarycount += 1
        #if stationarycount > 10:
        #    # if we're stuck, restart
        #    debugstr = 'HAD TO RESET! HAD TO RESET! HAD TO RESET!'
        #    stationarycount = 0
        #    for row in p:
        #        row = p0 + spread * np.random.rand(d)
        #lastmin = fmin
        frange = abs(np.max(fs) - np.min(fs))
        if isnan(frange) or isnan(fmin):
            raise RuntimeError('Downhill simplex donked out! calculated range is NaN.')
        debugstr += ' min value: ' + str(np.min(fs))
        debugstr += ' Simplex value range: '+str(frange)
        if tryno % 100 == 0:
            logging.debug(debugstr)
        if frange < rtol:
            satisfied = True
    bestind = np.argmin(fs)
    bestpt = p[bestind]
    bestmin = fs[bestind]
    return (bestpt, bestmin)

def rosenbrock(x):
    return (1 - x[0])*(1 - x[0]) + 100 * (x[1] - x[0]) * (x[1] - x[0])

# turn off debug messages
#logging.disable(logging.DEBUG)

def model_fit_mestimate(x, y, yfun, dist, p0, sigma=1, spread=1):
    """
    Computes a general M-estimate (maximum likelihood estimate)
    of a function yfun with parameters p to agiven set of data (x, y),
    whose y values have standard deviation sigma, in the presence of
    an arbitrary (normalized) error distribution

    en(x) = exp(-dist(x)).

    This is in contrast to a least-squares estimate, which assumes that
    the error is Gaussian (in other words, when mlfun = 0.5*x*x). This
    allows for a more robust fit of the model to the data: a Lorentzian,
    for example, will be more tolerant of outliers.

    :param x:
    :param y:
    :param fun:
    :param mfun:
    :param p0:
    :return:
    """

    # definition of the maximum likelihood function
    def mlfun(p):
        sm = np.sum(dist((y - yfun(x, p)) / sigma))
        #for k, (xi, yi) in enumerate(zip(x, y)):
        #    sm += dist((yi - yfun(xi, p)) / sigma)
        return sm

    # We'll use downhill simplex to minimize the maximum-likelihood
    # function. This is a somewhat slow method, but it's stable, and
    # guaranteed to converge to *something* eventually.
    # This line uses my homebrew downhill simplex. It sort of works.
    #pbest, bestmin = fmin_simplex(mlfun, p0, spread=spread)
    # This line uses Scipy's downhill simplex. It likely works better.
    optresult = scipy.optimize.minimize(mlfun, p0, method='Nelder-Mead')
    pbest = optresult.x
    # Additionally, we NEED to have a goodness-of-fit value to return,
    # but I don't know how to calculate that for M-estimate. So it's zero
    # for now.
    return pbest, 0.0

def rednoise(n, r=0.85):
    w = np.random.normal(size=n)
    redn = [w[0]]
    for j in range(len(w)-1):
        redj1 = r * redn[j] + np.sqrt(1 - r*r) * w[j+1]
        redn.append(redj1)
    return np.array(redn)

def dsin(x, p):
    return p[0] * np.sin(x*2*pi*p[1] + p[2])*np.exp(-x/p[3]) + p[4] + p[5]*np.exp(-x/p[6])

def dsinplus(x, p):
    return p[0] * np.cos(x*2*pi*p[1]*1E9)*np.exp(-x/p[2]*1E9) + p[3]

# testing our mestimate

#ydata = dsin(t, [1.5, 0.5, 1.0, 1.0, 0.75, 0.2, 2]) + 0.1*rednoise(len(t), r=0.85)

def ml_lorentz(z):
    return np.log(1 + 0.5 * z * z)

#bp, gof = model_fit_mestimate(t, ydata, dsin, ml_lorentz, [1, 1, 1, 1, 1, 1, 1], spread=0.3)
#yfit = dsin(t, bp)
#err = ydata - yfit

#print(bp)


#plt.plot(t, ydata)
#plt.plot(t, yfit)
#plt.hist(err, 50)
#plt.plot(t, err)
#plt.show()