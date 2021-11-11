'''
Created on 07.11.2021

@author: grustentier
'''

import os


# Returns the paths of all files that are located in a directory 
def collectFilePaths(directory):
    filePaths = []
    for root, _, files in os.walk(directory):  
        if root.endswith(os.sep) is False:root += os.sep
        filePaths.extend([root + file for file in files])
    return filePaths
