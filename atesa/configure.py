"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import mdtraj       # to support mdtraj calls in the input file
import numpy        # to support numpy  calls in input file
import numpy as np  # to support numpy  calls even if called as np
import sys
import os
import shutil
from jinja2 import Environment, FileSystemLoader
import typing
import pydantic

def configure(input_file, user_working_directory=''):
    """
    Configure the settings namespace based on the config file.

    Parameters
    ----------
    input_file : str
        Name of the configuration file to read
    user_working_directory : str
        User override for working directory (overrides value in input_file), ignored if set to ''

    Returns
    -------
    settings : argparse.Namespace
        Settings namespace object

    """

    class Settings(pydantic.BaseModel):
        # This class initializes the settings object with type hints. After being built, it gets exported as an
        # argparse.Namelist object, just for convenience.

        # Core settings required for all jobs and lacking defaults
        job_type: str
        batch_system: str
        restart: bool
        md_engine: str
        task_manager: str
        topology: str
        working_directory: str
        overwrite: bool

        # Batch template settings
        init_nodes: int = 1
        init_ppn: int = 1
        init_mem: str = '4000mb'
        init_walltime: str = '00:30:00'
        init_solver: str = 'sander'
        init_extra: str = ''
        prod_nodes: int = 1
        prod_ppn: int = 8
        prod_mem: str = '4000mb'
        prod_walltime: str = '02:00:00'
        prod_solver: str = 'sander'
        prod_extra: str = ''

        # File path settings (required for all jobs, but do have sensible defaults)
        path_to_input_files: str = os.path.dirname(os.path.realpath(__file__)) + '/data/input_files'
        path_to_templates: str = os.path.dirname(os.path.realpath(__file__)) + '/data/templates'

        # Required only for aimless shooting, equilibrium path sampling, and committor analysis
        cvs: typing.List[str] = ['']
        include_qdot: bool = True

        # Required only for aimless shooting and equilibrium path sampling
        initial_coordinates: typing.List[str] = ['']

        # Required only for aimless shooting and committor analysis
        commit_fwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])
        commit_bwd: typing.Tuple[typing.List[int], typing.List[int], typing.List[float], typing.List[str]] = ([-1], [-1], [-1], ['unset'])

        # Required only for committor analysis and equilibrium path sampling
        rc_definition: str = ''
        as_out_file: str = 'as_raw.out'
        rc_reduced_cvs: bool = True

        # Required only for aimless shooting
        min_dt: int = 1
        max_dt: int = 10
        always_new: bool = True
        resample: bool = False
        degeneracy: int = 1
        cleanup: bool = True
        information_error_checking: bool = True
        information_error_threshold: float = 0.1
        information_error_freq: int = 250
        information_error_override: bool = False
        information_error_max_dims: int = 6
        max_moves: int = -1     # also used by find_ts
        max_consecutive_fails: int = -1

        # Required only for committor analysis
        committor_analysis_n: int = 10
        committor_analysis_use_rc_out: bool = False
        path_to_rc_out: str = sys.path[0] + '/atesa/tests/test_data/rc.out'
        rc_threshold: float = 0.05

        # Required only for equilibrium path sampling
        eps_rc_min: float = -12
        eps_rc_max: float = 12
        eps_rc_step: float = 1
        eps_rc_overlap: float = 0.1
        eps_n_steps: int = 6
        eps_out_freq: int = 1
        eps_dynamic_seed: typing.Union[int, list] = 20  # int or list (int -> [int for window in eps_windows]; 0 or empty list turns off)
        samples_per_window: int = -1

        # Required only if restart = True
        restart_terminated_threads: bool = False

        # Not expected to be set by user
        DEBUG: bool = False     # True causes some functions to return dummy values for testing purposes
        pid: int = -1           # Process ID of information_error call used in aimless shooting
        information_error_overdue: bool = False    # Used for handling information_error calls cleanly
        dont_dump: bool = False     # When True, prevents dumping settings to settings.pkl

    # Import config file line-by-line using exec()
    try:
        lines = open(input_file, 'r').readlines()
    except FileNotFoundError:
        try:
            lines = open('atesa/' + input_file, 'r').readlines()     # for testing
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
    config_dict = {}
    config_dict.update(locals())
    settings = argparse.Namespace()
    settings.__dict__.update(Settings(**config_dict))

    # Override working directory if provided with user_working_directory
    if user_working_directory:
        settings.working_directory = user_working_directory

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

if __name__ == "__main__":
    configure('','')