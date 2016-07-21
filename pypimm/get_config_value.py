import re
import configparser

__author__ = 'alex'


def get_config_value(analysis, fp):
    cp = configparser.ConfigParser(inline_comment_prefixes='#;')
    with open(fp) as f:
        cp.read_file(f)
    analysis.new_configs(cp)
    return None