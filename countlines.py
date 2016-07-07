__author__ = 'alex'

import re, os
"""
This file counts the number of lines in all the python files in this directory.

2016-06-15: 1072.
2016-06-19: 1335
2016-06-24: 1716
"""

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

pythonFile = re.compile(r'.*\.py')

filesInDir = os.listdir('.')

pFilesInDir = []
for file in filesInDir:
    if pythonFile.match(file):
        pFilesInDir.append(file)

linecnt = 0

for file in pFilesInDir:
    linecnt += file_len(file)

print(pFilesInDir)

print(linecnt)