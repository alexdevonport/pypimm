import re
from PimmAnalysis import PimmAnalysis

__author__ = 'alex'


def get_config_value(analysis, fp):
    '''

    :param key: name of the config value we're looking up
    :param fp: filepath of the config file
    :return:
    '''
    f = open(fp, 'r')
    configstring = f.read()

    # the value corresponding to our key will either be a string or a number.
    #
    configre    = re.compile('([a-zA-Z ]*)=(\d+.?\d*)')

    for match in configre.findall(configstring):
        strkey = match[0]
        try:
            strval = float(match[1])
        except ValueError:
            strval = match[1]
        analysis.set_configs(strkey, strval)
    # That should be all the configs
    return None
