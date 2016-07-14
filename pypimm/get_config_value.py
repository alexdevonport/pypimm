import re
import configparser

__author__ = 'alex'


def get_config_value(analysis, fp):
    cp = configparser.ConfigParser(inline_comment_prefixes='#;')
    cp.read_file(open(fp))
    analysis.set_configs(cp)
    return None