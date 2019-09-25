"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import numpy        # to support numpy  calls in input file
import sys
import os
from jinja2 import Environment, FileSystemLoader

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

    # Import config file line-by-line using exec()
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

    # Set some default values
    if not settings.__contains__('DEBUG'):
        settings.DEBUG = False
    if not settings.__contains__('path_to_input_files'):
        settings.path_to_input_files = sys.path[0] + '/data/input_files'
    if not settings.__contains__('path_to_templates'):
        settings.path_to_templates = sys.path[0] + '/data/templates'

    # Format directories properly (no trailing '/')
    if settings.working_directory[-1] == '/':
        settings.working_directory = settings.working_directory[:-1]
    if settings.path_to_input_files[-1] == '/':
        settings.path_to_input_files = settings.path_to_input_files[:-1]
    if settings.path_to_templates[-1] == '/':
        settings.path_to_templates = settings.path_to_templates[:-1]
    try:    # this setting is optional
        if settings.path_to_rc_out[-1] == '/':
            settings.path_to_rc_out = settings.path_to_rc_out[:-1]
    except AttributeError:
        pass

    # Set Jinja2 environment
    if os.path.exists(settings.path_to_templates):
        settings.env = Environment(loader=FileSystemLoader(settings.path_to_templates))
    else:
        sys.exit('Error: could not locate templates folder: ' + settings.path_to_templates)

    # Check that each given value is valid # todo: implement this for every option
    if not settings.job_type in ['aimless_shooting', 'committor_analysis', 'equilibrium_path_sampling', 'isee']:
        raise ValueError('unsupported job_type: ' + str(settings.job_type))
    if not type(settings.initial_coordinates) == list:
        raise ValueError('initial_coordinates must be provided as a list (even if it is of length 1)')
    for item in settings.initial_coordinates:
        if '/' in item:
            raise ValueError('file names in initial_coordinates should not contain \'/\' characters (they should exist in the same directory from which ATESA is called)')


    return settings
