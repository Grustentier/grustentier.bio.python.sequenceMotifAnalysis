'''
Created on 07.11.2021

@author: grustentier
'''

import re


def cleanString(string, replacement='-'):
    a = re.sub('[^a-zA-Z0-9.?]', replacement, string) 
    return re.sub(replacement + '+', replacement, a)


def asBoolean(s):
    if str(s).lower() not in ['false', 'true', '1', '0']:
        raise ValueError('Not a valid boolean string')
    return str(s).lower() == 'true' or str(s).lower() == '1'
