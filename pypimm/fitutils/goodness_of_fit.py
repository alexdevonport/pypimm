import numpy as np
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

def redchi2(fun, x, y, p, sigma=1.0):
    n = np.size(y)
    m = np.size (p)
    return np.sum(np.power((y - fun(x, *p)) / (sigma * (n - m)),2))

