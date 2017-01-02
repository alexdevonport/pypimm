__author__ = 'alex'
import json

class PimmAnalysis():
    def __init__(self, name=None, timebase=None, raw_data=None, metadata=None,
                 fits=None, results=None, progress=None, configs=None):
        self.name = name
        self.timebase = timebase
        self.raw_data = raw_data
        self.metadata = metadata
        self.fits = fits
        self.rawfits = None
        self.results = results
        self.progress = progress
        self.configs = configs
        self.fields = None

    def set_name(self, newname):
        self.name = newname
        return None

    def get_name(self):
        return self.name

    def set_timebase(self, newtimebase):
        self.timebase = newtimebase
        return None

    def get_timebase(self):
        return self.timebase

    def set_fields(self, newfields):
        self.fields = newfields

    def get_fields(self):
        return self.fields

    def set_raw_data(self, newrawdata):
        self.raw_data = newrawdata
        return None

    def get_raw_data(self):
        return self.raw_data

    def set_metadata(self, newmetedata):
        self.metadata = newmetedata
        return None

    def get_metadata(self):
        return self.metadata

    def print_metadata(self):
        print(json.dumps(self.metadata, indent=4))

    def set_fits(self, newfits):
        self.fits = newfits
        return None

    def get_fits(self):
        return self.fits

    def set_rawfits(self, newfits):
        self.rawfits = newfits
        return None

    def get_rawfits(self):
        return self.rawfits


    def set_results(self, newresults):
        self.results = newresults
        return None

    def get_results(self):
        return self.results

    def set_progress(self, newprogress):
        self.progress = newprogress
        return None

    def get_progress(self):
        return self.progress

    def new_configs(self, newconfigs):
        self.configs = newconfigs
        return None

    def set_configs(self, section, option, value):
        self.configs.set(section, option, value)
        return None


    def get_configs(self, configkey=None):
        if configkey is None:
            return self.configs
        else:
            ckey1, ckey2 = configkey
            return self.configs[ckey1][ckey2]