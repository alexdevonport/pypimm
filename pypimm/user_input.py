import re
__author__ = 'alex'

def user_input(analysis):
    print('Go and look at the fit waveforms now.')
    print('Are there any you don\'t want in the fit?')
    fits = analysis.get_fits()
    fitkeys = list(analysis.get_fits().keys())

    for index, key in enumerate(fitkeys):
        print(index, key)

    print('Select which ones you want to ignore as a list of indices,')
    print('e.g. \'0, 1, \' to remove {} and {}'.format(fitkeys[0], fitkeys[1]))
    s = input('> ')

    for delind in str_to_ints(s):
        del fits[fitkeys[delind]]

    print('Good choice! We will only fit the following for Ms and Hk.')
    print(fits.keys())

    return None

def str_to_ints(s):
    s = re.sub('[\[\]\(\)\,]', ' ', s)
    return [int(elt) for elt in s.split()]
