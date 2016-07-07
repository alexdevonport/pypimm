__author__ = 'alex'
import shutil, time


class ProgBar():
    def __init__(self, msg=None, maxn=None, w=None):
        self.msg = msg
        self.maxn = maxn
        self.n = 0
        self.w = w
        self.update(0)

    def set_msg(self, newmsg):
        self.msg = newmsg
        return None

    def set_maxn(self, newmaxn):
        self.maxn = newmaxn
        return None

    def set_w(self, neww):
        self.w = neww
        return None

    def update(self, newn=None):
        if newn != None:
            self.n += newn
        nhash = int(self.n / self.maxn * (self.w - 2))
        # Generate a zero-padded message string
        cols = shutil.get_terminal_size((80, 24))[0]
        lmsg = cols - self.w - 3

        # msg will contain all test before the progress bar, including the white
        # space. We need to take care of two cases: truncating the message,
        # or padding it with spaces.
        if len(self.msg) <= lmsg:
            msg = self.msg + ' ' * (lmsg - len(self.msg))
        else:
            msg = self.msg[:lmsg-3] + '...'
        pbar = '[' + '#' * nhash + '-' * (self.w - nhash - 2) + ']'
        pbstr = '\r' + msg + '   ' + pbar
        print(pbstr, end='')

    def restart(self, newmsg=None):
        if newmsg != None:
            self.msg = newmsg
        self.n = 0
        print('\n')  # Start new progress bar on next line
        self.update(0)

