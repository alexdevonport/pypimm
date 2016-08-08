__author__ = 'alex'
import scipy.stats
import matplotlib.pyplot as plt
import os

def error_analysis(timebase, signal, fit, name):
    """
    This function is used to take a look at the deviation between PyPIMM's
    best fit for the signal and the signal itself. It does the following:
        -compute the histogram and spectrum of the fit error signal (data - fit)
        -attempt to fit the histogram to three error distributions (Gaussian
         single-sided exponential, Lorentzian)
        -Attempts to fit the noise spectrum to several noise power spectral densities
    :param timebase:
    :param signal:
    :param analysis:
    :return:
    """

    errsig = signal - fit

    plt.clf()
    plt.plot(timebase*1E9, errsig)
    plt.xlabel('time (ns)')
    plt.ylabel('signal - best fit (V)')
    plt.title('Best Fit Residue for '+name)
    plt.savefig(os.path.join('.','Error','time-domain','residue-'+name+'-time.png'))
    plt.clf()
    scipy.stats.probplot(errsig, plot=plt)
    plt.title('Best-fit residue histogram for '+name)
    plt.savefig(os.path.join('.','Error','normal-plot','residue-'+name+'-normal.png'))
    plt.clf()

    return None