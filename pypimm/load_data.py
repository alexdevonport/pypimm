import csv
import re
import pandas as pd
import json
from collections import OrderedDict
import numpy as np

def load_data(a, fp):
    if fp.endswith('.json'):
        t, sigs, md = load_data_json(fp)
        a.set_metadata(md)
        a.print_metadata()
    elif fp.endswith('.txt'):
        t, sigs = load_data_csv(fp)
    else:
        raise Exception('Data file format not recognized. JSON or TXT.')
    a.set_timebase(t)
    a.set_raw_data(sigs)

def load_data_csv(filepath):
    timebase_key = 'time(s)'
    r = pd.read_csv(filepath, sep='\t')
    signal_header = re.compile(r'.*signal.*')
    signals = OrderedDict()
    for key in r.keys():
        if signal_header.match(key):
            key_nospace = key.replace(' ', '-')
            signals[key_nospace] = r[key]
    return r[timebase_key], signals

def load_data_json(filepath):
    with open(filepath) as f:
        signals = OrderedDict()
        j = json.load(f)
        tb = []
        metadata = j[' metadata']
        for pdata in j['raw data']:
            tb = np.array(pdata['timebase'])
            sname = '{:f}-Oe-signal'.format(pdata['bias field'])
            signals[sname] = np.array(pdata['signal'])
        return tb, signals, metadata
