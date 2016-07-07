__author__ = 'alex'

class PimmAnalysis():
    def __init__(self, name=None, timebase=None, raw_data=None,
                 fits=None, results=None, progress=None, configs=None):
        self.name = name
        self.timebase = timebase
        self.raw_data = raw_data
        self.fits = fits
        self.results = results
        self.progress = progress
        self.configs = configs

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

    def set_raw_data(self, newrawdata):
        self.raw_data = newrawdata
        return None

    def get_raw_data(self):
        return self.raw_data

    def set_fits(self, newfits):
        self.fits = newfits
        return None

    def get_fits(self):
        return self.fits

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

    def set_configs(self, newconfigkey = None, newconfigval=None):
        if self.configs is None:
            self.configs = {newconfigkey:newconfigval}
        else:
            self.configs[newconfigkey] = newconfigval
        return None

    def get_configs(self, configkey=None):
        if configkey is None:
            return self.configs
        else:
            return self.configs[configkey]