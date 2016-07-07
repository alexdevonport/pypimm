__author__ = 'alex'


def pimm_write_xls(name, data, fits, properties):
    """
    Creates an Excel file containing the results of the analysis that's just
    been done. It has three sheets:

    (1) the raw PIMM signal data
    (2) signal spectra and best-fit curves
    (3) final results (ms, hk, damping, etc.

    :param name: what to name the data file
    :param data: raw PIMM data from file
    :param fits: fit parameters gained from
    :param properties: material properties found from fits
    :return:
    """
    return None