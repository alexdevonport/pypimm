import csv
import re
import pandas as pd
import json
from collections import OrderedDict
import numpy as np
import pprint

def load_data(a, fp):
    if fp.endswith('.json'):
        t, sigs, hs, md = load_data_json(fp)
        a.set_metadata(md)
        a.print_metadata()
    elif fp.endswith('.txt'):
        t, sigs, hs = load_data_csv(fp)
    else:
        raise Exception('Data file format not recognized. JSON or TXT.')
    a.set_timebase(t)
    a.set_raw_data(sigs)
    a.set_fields(hs)

def load_data_csv(filepath):
    timebase_key = 'time(s)'
    hs = []
    signames = []
    r = pd.read_csv(filepath, sep='\t')
    signal_header = re.compile(r'([+-]?\d+\.?\d*) Oe signal')
    signals = OrderedDict()
    # get all signals (not biases or sats) from the csv
    for key in r.keys():
        mo = signal_header.match(key)
        if not mo == None:
        #if signal_header.match(key):
            key_nospace = key.replace(' ', '-')
            hs.append(float(mo.group(1)))
            signames.append(key_nospace)
            signals[key_nospace] = r[key]
    # now sort them
    hs_s, signames_s = sortwith(hs, signames)
    signals_s = OrderedDict()
    for key in signames:
        signals_s[key] = signals[key]
    return r[timebase_key], signals_s, hs_s

def load_data_json(filepath):
    with open(filepath) as f:
        signals = OrderedDict()
        hs = []
        j = json.load(f)
        tb = []
        metadata = j[' metadata']
        for pdata in j['raw data']:
            hs.append(pdata['bias field'])
            tb = np.array(pdata['timebase'])
            sname = '{:f}-Oe-signal'.format(pdata['bias field'])
            signals[sname] = np.array(pdata['signal'])
        return tb, signals, hs, metadata

def sortwith(a, *b):
    """
    Sorts the sequence a, and performs the same rearrangements
    on all the sequences in the tuple b. That is, the sequences
    in b are sorted "along with" a.
    :param a: sequence
    :param b: tuple of sequences
    :return: tuple of sequences. Returns all of the sorted
    sequences, starting with a, then all the ones in b.
    """
    r = (list(t) for t in zip(*sorted(zip(a, *b))))
    return r