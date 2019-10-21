"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import numpy        # to support numpy  calls in input file
import numpy as np  # to support numpy  calls even if called as np
import sys
import os
import shutil
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
    for line in lines:      # each line in the input file is just python code setting a variable;
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

    # Check that each given value is valid
    if not type(settings.job_type) == str:
        raise ValueError('job_type must be given as a string')
    if not type(settings.batch_system) == str:
        raise ValueError('batch_system must be given as a string')
    if not type(settings.restart) == bool:
        raise ValueError('restart must be given as a boolean')
    if not type(settings.md_engine) == str:
        raise ValueError('md_engine must be given as a string')
    if not type(settings.task_manager) == str:
        raise ValueError('task_manager must be given as a string')
    if not type(settings.topology) == str:
        raise ValueError('topology must be given as a string')
    if not type(settings.working_directory) == str:
        raise ValueError('working_directory must be given as a string')
    if not type(settings.overwrite) == bool:
        raise ValueError('overwrite must be given as a boolean')
    if not type(settings.path_to_templates) == str:
        raise ValueError('path_to_templates must be given as a string')
    if not type(settings.path_to_input_files) == str:
        raise ValueError('path_to_input_files must be given as a string')
    if settings.job_type in ['aimless_shooting', 'equilibrium_path_sampling']:
        if not type(settings.init_nodes) == int:
            raise ValueError('init_nodes must be given as an integer')
        if not type(settings.init_ppn) == int:
            raise ValueError('init_ppn must be given as an integer')
        if not type(settings.init_mem) in [str, int]:
            raise ValueError('init_mem must be given as a string or integer')
        if not type(settings.init_walltime) == str:
            raise ValueError('init_walltime must be given as a string')
        if not type(settings.init_solver) == str:
            raise ValueError('init_solver must be given as a string')
        if not type(settings.initial_coordinates) == list:
            raise ValueError('initial_coordinates must be given as a list (even if it is of length 1)')
        if not type(settings.commit_fwd) == list:
            raise ValueError('commit_fwd must be given as a list (even if it is of length 1)')
        for item in settings.commit_fwd:
            if not type(item) == list:
                raise ValueError('items in commit_fwd must be given as lists (even if they are of length 1)')
        for item in settings.commit_fwd[0]:
            if not type(item) == int:
                raise ValueError('items in the first list within commit_fwd must be given as integers')
        for item in settings.commit_fwd[1]:
            if not type(item) == int:
                raise ValueError('items in the second list within commit_fwd must be given as integers')
        for item in settings.commit_fwd[2]:
            if not type(item) in [int, float]:
                raise ValueError('items in the third list within commit_fwd must be given as integers or floats')
        for item in settings.commit_fwd[3]:
            if not item in ['gt', 'lt']:
                raise ValueError('items in the fourth list within commit_fwd must be either \'gt\' or \'lt\'')
        if not type(settings.commit_bwd) == list:
            raise ValueError('commit_bwd must be given as a list (even if it is of length 1)')
        for item in settings.commit_bwd:
            if not type(item) == list:
                raise ValueError('items in commit_bwd must be given as lists (even if they are of length 1)')
        for item in settings.commit_bwd[0]:
            if not type(item) == int:
                raise ValueError('items in the first list within commit_bwd must be given as integers')
        for item in settings.commit_bwd[1]:
            if not type(item) == int:
                raise ValueError('items in the second list within commit_bwd must be given as integers')
        for item in settings.commit_bwd[2]:
            if not type(item) in [int, float]:
                raise ValueError('items in the third list within commit_bwd must be given as integers or floats')
        for item in settings.commit_bwd[3]:
            if not item in ['gt', 'lt']:
                raise ValueError('items in the fourth list within commit_bwd must be either \'gt\' or \'lt\'')
    if settings.job_type in ['aimless_shooting', 'equilibrium_path_sampling', 'committor_analysis']:
        if not type(settings.prod_nodes) == int:
            raise ValueError('prod_nodes must be given as an integer')
        if not type(settings.prod_ppn) == int:
            raise ValueError('prod_ppn must be given as an integer')
        if not type(settings.prod_mem) in [str, int]:
            raise ValueError('prod_mem must be given as a string or integer')
        if not type(settings.prod_walltime) == str:
            raise ValueError('prod_walltime must be given as a string')
        if not type(settings.prod_solver) == str:
            raise ValueError('prod_solver must be given as a string')
        if not type(settings.cvs) == list:
            raise ValueError('cvs must be given as a list (even if it is of length 1)')
        for item in settings.cvs:
            if not type(item) == str:
                raise ValueError('items in the cvs list must be given as strings')
        if not type(settings.include_qdot) == bool:
            raise ValueError('include_qdot must be given as a boolean')
    if settings.job_type in ['equilibrium_path_sampling', 'committor_analysis']:
        if not type(settings.rc_definition) == str:
            raise ValueError('rc_definition must be given as a string')
        if not type(settings.as_out_file) == str:
            raise ValueError('as_out_file must be given as a string')
        if not type(settings.rc_reduced_cvs) == bool:
            raise ValueError('rc_reduced_cvs must be given as a boolean')
    if settings.job_type == 'aimless_shooting':
        if not type(settings.min_dt) == int:
            raise ValueError('min_dt must be given as an integer')
        if not type(settings.max_dt) == int:
            raise ValueError('max_dt must be given as an integer')
        if not type(settings.always_new) == bool:
            raise ValueError('always_new must be given as a boolean')
        if not type(settings.resample) == bool:
            raise ValueError('resample must be given as a boolean')
        if not type(settings.degeneracy) == int:
            raise ValueError('degeneracy must be given as an integer')
        if not type(settings.information_error_checking) == bool:
            raise ValueError('information_error_checking must be given as a boolean')
        if settings.information_error_checking:
            if not type(settings.information_error_freq) == int:
                raise ValueError('information_error_freq must be given as an integer')
            if not type(settings.information_error_override) == bool:
                raise ValueError('information_error_override must be given as a boolean')
            if not settings.information_error_override and len(settings.cvs) < 5:
                raise RuntimeError('information_error_checking = True requires that the cvs option has at least five entries to'
                                   ' be able to dynamically identify suitable reaction coordinates for evaluation of the model '
                                   'information error, and many more than five is strongly recommended for higher confidence in'
                                   ' the results. To override this behavior (only to be used if you are extremely confident '
                                   'that your set of CVs is complete), set information_error_override = True')
    if settings.job_type == 'committor_analysis':
        if not type(settings.committor_analysis_n) == int:
            raise ValueError('committor_analysis_n must be given as an integer')
        if not type(settings.committor_analysis_use_rc_out) == bool:
            raise ValueError('committor_analysis_use_rc_out must be given as a boolean')
        if settings.committor_analysis_use_rc_out:
            if not type(settings.path_to_rc_out) == str:
                raise ValueError('path_to_rc_out must be given as a string')
            if not os.path.exists(settings.path_to_rc_out):
                raise RuntimeError('reaction coordinate output (rc_out) file not found: ' + settings.path_to_rc_out)
            if not type(settings.rc_threshold) in [int, float]:
                raise ValueError('rc_threshold must be given as an integer or float')
    if settings.job_type in ['aimless_shooting', 'equilibrium_path_sampling'] or (settings.job_type == 'committor_analysis' and not settings.committor_analysis_use_rc_out):
        filenames = []  # check for duplicates in initial_coordinates
        for item in settings.initial_coordinates:
            if not type(item) == str:
                raise ValueError('entries in the initial_coordinates list must be given as strings')
            if '/' in item:
                filenames.append(item[item.rindex('/') + 1:])
            else:
                filenames.append(item)
        seen = []
        for filename in filenames:
            if filename in seen:
                raise RuntimeWarning(
                    'duplicate filename in initial_coordinates: ' + filename + '\nOnly the right-most file'
                                                                               ' with this name (independent of path) will be used to seed each of the threads whose '
                                                                               'initial coordinates share this name')
            seen.append(filename)
    if settings.job_type == 'equilibrium_path_sampling':
        if not type(settings.eps_rc_min) in [int, float]:
            raise ValueError('eps_rc_min must be given as an integer or float')
        if not type(settings.eps_rc_max) in [int, float]:
            raise ValueError('eps_rc_max must be given as an integer or float')
        if not type(settings.eps_rc_step) in [int, float]:
            raise ValueError('eps_rc_step must be given as an integer or float')
        if not type(settings.eps_rc_overlap) in [int, float]:
            raise ValueError('eps_rc_overlap must be given as an integer or float')
        if not type(settings.eps_n_steps) == int:
            raise ValueError('eps_n_steps must be given as an integer')
        if not type(settings.eps_out_freq) == int:
            raise ValueError('eps_out_freq must be given as an integer')
        if not type(settings.eps_dynamic_seed) in [int, list]:
            raise ValueError('eps_dynamic_seed must be given as an integer or list')
    if settings.restart:
        if not type(settings.restart_terminated_threads) == bool:
            raise ValueError('restart_terminated_threads must be given as a boolean')

    # Initialize EPS settings based on contents of config file
    if settings.job_type == 'equilibrium_path_sampling':
        eps_lower_boundaries = numpy.arange(settings.eps_rc_min, settings.eps_rc_max, settings.eps_rc_step)
        settings.eps_bounds = []
        for lower_bound in eps_lower_boundaries:
            settings.eps_bounds.append([lower_bound - settings.eps_rc_overlap, lower_bound + settings.eps_rc_step + settings.eps_rc_overlap])

        if settings.eps_dynamic_seed:
            if isinstance(settings.eps_dynamic_seed, int):
                settings.eps_dynamic_seed = [settings.eps_dynamic_seed for null in settings.eps_bounds]
            settings.eps_empty_windows = settings.eps_dynamic_seed

    return settings
