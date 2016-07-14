import numpy as np
import itertools as it
from math import pi
import matplotlib.pyplot as plt
__author__ = 'alex'


def grid_lsq(fun, x, y, pdims):
    bestp = []
    bestres = 1E100
    def res(rfun, rp):
        plt.plot(x, rfun(x, *rp))
        #plt.show()
        return np.sum(np.power(y - rfun(x, *rp), 2))
    plists = [np.linspace(p0, p1, n) for (p0, p1, n) in pdims]
    pspace = it.product(*plists)
    for p in pspace:
        pres = res(fun, p)
        if pres < bestres:
            bestres = pres
            bestp = p
    return bestp

def main():

    def dsin(x, p0, p1, p2):
        return p0 * np.sin(2 * pi * p1 * x) * np.exp(-x * p2)
    x = np.linspace(0, 10, 1000)
    y = dsin(x, 1, 2, 3)
    pguess = grid_lsq(dsin, x, y, [(0.1, 2, 10), (1, 5, 10), (0, 5, 10)])
    print(pguess)
    plt.clf()
    plt.plot(x, y, x, dsin(x, *pguess))
    plt.show()
    print(pguess)

if __name__ == '__main__':
    main()