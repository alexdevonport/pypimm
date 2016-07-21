#!/usr/bin/env python

import os
import sys
import logging

from PimmAnalysis import  PimmAnalysis
from load_data import load_data
from fitutils import fit_data, characterize_fits
from generate_report import generate_report
from get_config_value import get_config_value
from user_input import user_input

#logging.basicConfig(level=logging.DEBUG)
logging.disable(logging.DEBUG)
logging.disable(logging.WARNING)
logging.basicConfig(format=' %(asctime)s - %(levelname)s - %(message)s',
                    filename='skipped_data.log',
                    level=logging.WARNING,)

timebase_key = 'time(s)'

# initialize data structures for signal data, signal fits, and characterization
data = {} # contains pimm signals & timebase
fits = {} # contains fit parameters (freq, damping, start point) for each signal
properties = {} # contains final analysis results (Ms, damping sweep, Hk, rm, etc)

# load PIMM data file from given file path
# accepts multiple data files, and do PIMM analysis
# on each one and put the results in a different folder

# check to see if we got a filename at all
if len(sys.argv) < 2:
    # TODO: replace warning message & exit with a file dialog
    print('Did you supply a file?')
    sys.exit()

pypimmdir = os.path.dirname(os.path.abspath(__file__))

#parse options
if '--handpick' in sys.argv:
    handpick = 'yes'
else:
    handpick = 'no'

for filepath in sys.argv[1:]:
    analysis = PimmAnalysis()
    fpnoext = os.path.split(filepath)[1][:-4] # takes the path and returns just the file name
    fpnoext = fpnoext.replace(' ','')
    fpnoext = fpnoext.replace('_','-')
    analysis.set_name(fpnoext)

    # Get the signals ant timebase out of the data file and into the object
    try:
        load_data(analysis, filepath)
    except OSError:
        continue
    # Get configuration values from the config file
    configfp = os.path.join(pypimmdir, 'pypimmconfig.txt')
    get_config_value(analysis, configfp)
    analysis.set_configs('fit-params', 'hand pick values', handpick)

    # Make directories to store results in. It may be a good idea to
    # replace this with a function that accepts a 'directory structure'
    # member of the Analysis object and generates all the directories from
    # that, some time down the road.
    # make a directory named after the data file and cd into it
    os.mkdir(os.path.join('.', fpnoext))
    os.chdir(os.path.join('.', fpnoext))
    # in the new data directory, make data folders
    os.mkdir(os.path.join('.', 'spectra'))
    os.mkdir(os.path.join('.', 'envelopes'))
    os.mkdir(os.path.join('.', 'sigfits'))
    os.mkdir(os.path.join('.', 'report'))
    os.mkdir(os.path.join('.', 'Error'))
    os.mkdir(os.path.join('.', 'Error','time-domain'))
    os.mkdir(os.path.join('.', 'Error','histogram'))

    # for each signal, calculate fit parameters
    fit_data(analysis)

    if handpick == 'yes':
        user_input(analysis)

    # with all signal fits done, calculate material properties
    characterize_fits(analysis)
    # TODO: compile results from each previous stage into a report file
    generate_report(analysis, pypimmdir)

    # TODO: now that the analysis for this file is over, cd to the original directory
    os.chdir('..')

# TODO: Now that ALL analysis is over, create a 'highlights' document showing
# all frequency and damping plots

# Now that everything is done, report to the user and quit.
print('Analysis complete.')
