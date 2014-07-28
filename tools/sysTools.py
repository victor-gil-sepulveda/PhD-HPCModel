'''
Created on 18/02/2014

@author: victor
'''

import os
import errno
import json

def createFolder(directory_path, ensure_writability = False):
    """
    Creates a directory (with subdirs) if it doesn't exist.

    @param directory_path: the path of the directory and subdirectories to be created.

    @return: True if it was possible to create the folder, False otherwise.
    """
    if ensure_writability:
        if not os.access(os.path.dirname(directory_path), os.W_OK):
            return False
    try:
        os.makedirs(directory_path)
        return True
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    return False

def saveDicInJSONFile(filename, dic):
    """
    Saves a dictionary in a file in JSON format.

    @param filename: is the path with the filename where we want to store the data.
    """
    open(filename, "w").write(json.dumps(dic, sort_keys = False, indent = 4, separators = (',', ': ')))

def loadDicInJSONFile(filename):
    """
    Loads a dictionary from a JSON file.

    @param filename: is the path with the filename where we want to store the data.

    @return: A Dictionary with the data.
    """
    return json.loads(open(filename,"r").read())


