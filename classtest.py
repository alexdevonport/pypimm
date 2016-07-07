__author__ = 'alex'

class a():

    def __init__(self, av):
        self.av = av
        return None

    def set_av(self, newav):
        self.av = newav
        return None

    def get_av(self):
        return self.av

def afun(ai, v):
    ai.set_av(v)
    return None

ainst = a(3)
print(ainst.get_av())
afun(ainst, 100)
print(str(ainst.get_av()))