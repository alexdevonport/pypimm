import csv
import re
import pandas as pd
from collections import OrderedDict

def load_data(analysis, filepath):
    timebase_key = 'time(s)'
    r = pd.read_csv(filepath, sep='\t')
    signal_header = re.compile(r'.*signal.*')
    signals = OrderedDict()
    for key in r.keys():
        if signal_header.match(key):
            key_nospace = key.replace(' ', '-')
            signals[key_nospace] = r[key]
    # That's all the signals in one dict; we'll add that to the analyis object
    analysis.set_raw_data(signals)
    # add timebase as well
    analysis.set_timebase(r[timebase_key])
    return None

"""
def load_data(analysis, filepath):
    '''
    This function loads the data from a PIMM signal sweep into a data
    structure used by the rest of the program. The data file is plaintext
    delimited by tabs and newlines, so the csv library is perfect to use here.

    :rtype : dict
    :param filepath: location of the data folder.
    :return: dict containing data vectors, keyed as either 'timebase' for the
    time vector, or as 'x Oe' for the x Oe signal data. Additionally, it might be
    nice to later add additional entries with other data that could be gleaned from the
    data file, e.g. saturation field, number of averages, etc.
    '''

    # TODO: Stop using dicts, and use ordered dicts instead.
    timebase_key = 'time(s)'
    r = OrderedDict([])
    with open(filepath, 'r') as csvfile:
        for row in csv.DictReader(csvfile, delimiter='\t'):
            # iterate over evert row in the csv file. The data cell
            # in the present row under the column with header x
            # will be row['x'].
            for key, val in row.items():
                # for the first row, initialize the key with the first value
                if(key not in r.keys()):
                    try:
                        r[key] = [float(val)]
                    except TypeError:
                        r[key] = [0.0]
                # after that, we can just append to the list value for each key
                else:
                    try:
                        r[key].append(float(val))
                    except TypeError:
                        r[key].append(0.0)
    # Now that we have all the columns from the data file, we'll
    # select the ones we want to keep
    signal_header = re.compile(r'.*signal.*')
    signals = {}
    for key in r.keys():
        if signal_header.match(key):
            key_nospace = key.replace(' ', '-')
            signals[key_nospace] = r[key]
    # That's all the signals in one dict; we'll add that to the analyis object
    analysis.set_raw_data(signals)
    # add timebase as well
    analysis.set_timebase(r[timebase_key])

    return None
"""
