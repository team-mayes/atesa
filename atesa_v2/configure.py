"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import numpy        # to support numpy  calls in input file
import os

def configure(input_file):
    os.getcwd()
    try:
        lines = open(input_file, 'r').readlines()
    except FileNotFoundError:
        try:
            lines = open('atesa_v2/' + input_file, 'r').readlines()     # for testing
        except:
            lines = open(input_file, 'r').readlines()   # to reproduce original error
    for line in lines:  # each line in the input file is just python code setting a variable;
        exec(line)      # this means that comments are supported using '#' and whitespace is ignored.

    # Define settings namespace to store all these variables
    settings = argparse.Namespace()
    settings.__dict__.update(locals())

    return settings
