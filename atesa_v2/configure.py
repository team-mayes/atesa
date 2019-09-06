"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import numpy        # to support numpy  calls in input file
import sys

def configure(input_file):
    """
    Configure the settings namespace based on the config file.

    Parameters
    ----------
    input_file : str
        Name of the configuration file to read

    Returns
    -------
    settings : argparse.Namespace
        Settings namespace object

    """

    # Set some default values
    DEBUG = False
    path_to_input_files = sys.path[0] + '/atesa_v2/data/input_files'
    path_to_templates = sys.path[0] + '/atesa_v2/data/templates'

    try:
        lines = open(input_file, 'r').readlines()
    except FileNotFoundError:
        try:
            lines = open('atesa_v2/' + input_file, 'r').readlines()     # for testing
        except:
            lines = open(input_file, 'r').readlines()   # to reproduce original error
    line_index = 0
    for line in lines:  # each line in the input file is just python code setting a variable;
        line_index += 1
        try:
            exec(line)      # this means that comments are supported using '#' and whitespace is ignored.
        except Exception as e:
            raise ValueError('error raised while reading line ' + str(int(line_index)) + ' of configuration file '
                             + input_file + ': ' + str(e))

    # Define settings namespace to store all these variables
    settings = argparse.Namespace()
    settings.__dict__.update(locals())

    return settings
