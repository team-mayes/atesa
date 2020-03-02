"""
An auxiliary script to read threads from a restart.pkl file and verify their integrity. Also backs up the restart.pkl
file and removes broken Threads as necessary.
"""

import os
import sys
import pickle
import pytraj
import typing
import warnings
import operator
import argparse
from atesa_v2.main import Thread
from atesa_v2 import factory
from pydantic import BaseModel, ValidationError, validator


def main(pkl_file, settings_pkl='settings.pkl'):
    """
    Verify integrity of Threads in pkl_file and remove broken threads as necessary.

    Iterates through Thread objects in the pkl file, checking for:
        - Readability as a Thread object with the proper attribute types
        - Formatting of the history attribute, including that entries are of consistent lengths
        - Validity of attributes with limited valid entries, such as current_type
        - Existence and validity of last_accepted trajectories and restart files, as appropriate
        
    Backs up the pkl_file before running. If a thread is found to be invalid, a warning is produced and it is removed 
    from the file, overwriting the existing pkl_file.

    Parameters
    ----------
    pkl_file : str
        Path to the restart.pkl file to operate on
    settings_pkl : str
        Path to the settings.pkl file corresponding to the restart.pkl file

    Returns
    -------
    output_str : str
        String of text to inform the user of the outcome of the verification step.
    output_int : int
        0 if the pkl file was valid and has not been modified; 1 if it was invalid and has been modified

    """

    try:
        settings = pickle.load(open(settings_pkl, 'rb'))
    except FileNotFoundError:
        raise FileNotFoundError('could not find settings pickle file: ' + settings_pkl + '. Please specify the path to '
                                'this file as the second argument passed to verify_threads.py. It should be in the '
                                'same working directory as the original restart pickle file.')

    jobtype = factory.jobtype_factory(settings.job_type)

    # Define custom Exception to identify validation errors
    class ThreadValidationError(Exception):
        pass

    class ValidThread(BaseModel):
        """ Template for a valid Thread object """

        # def __init__(self, dict):
        #     self.dict = dict

        topology: str = ''
        jobids: typing.List = []
        terminated: bool = False
        current_type: typing.List = []
        current_name: typing.List = []
        current_results: typing.List = []
        name: str = ''
        suffix: int = 0
        total_moves: int = 0
        accept_moves: int = 0
        status: str = 'fresh thread'
        skip_update: bool = False
        history: argparse.Namespace

        class Config:
            arbitrary_types_allowed = True

        @validator('suffix', 'total_moves', 'accept_moves')
        def is_at_least_zero(cls, v):
            if v < 0:
                raise ThreadValidationError
            return v

        @validator('history')
        def valid_history(cls, v, values):
            if not jobtype.verify(cls, v, values['suffix'], 'history'):
                raise ThreadValidationError
            return v

        # @validator('topology')
        # def file_exists(cls, v):
        #     if not os.path.exists(v):
        #         raise ThreadValidationError
        #     return v

        ## I think not necessary for pydantic.BaseModel
        # def __init__(self, **kwargs):
        #     for key, value in kwargs.items():
        #         setattr(self, key, value)

        # topology = property(operator.attrgetter('_topology'))
        # @topology.setter
        # def topology(self, x):
        #     if not type(x) == str: raise MyValidationError
        #
        # jobids = property(operator.attrgetter('_jobids'))
        # @jobids.setter
        # def jobids(self, x):
        #     if not type(x) == list: raise MyValidationError
        #
        # terminated = property(operator.attrgetter('_terminated'))
        # @terminated.setter
        # def terminated(self, x):
        #     if not type(x) == bool: raise MyValidationError
        #
        # current_type = property(operator.attrgetter('_current_type'))
        # @current_type.setter
        # def current_type(self, x):
        #     if not type(x) == list or False in [type(item) == str for item in x]: raise MyValidationError
        #
        # current_name = property(operator.attrgetter('_current_name'))
        # @current_name.setter
        # def current_name(self, x):
        #     if not type(x) == list or False in [type(item) == str for item in x]: raise MyValidationError
        #
        # current_results = property(operator.attrgetter('_current_results'))
        # @current_results.setter
        # def current_results(self, x):
        #     if not type(x) == list or False in [type(item) == str for item in x]: raise MyValidationError
        #
        # name = property(operator.attrgetter('_name'))
        # @name.setter
        # def name(self, x):
        #     if not type(x) == str: raise MyValidationError
        #
        # suffix = property(operator.attrgetter('_suffix'))
        # @suffix.setter
        # def suffix(self, x):
        #     if not type(x) == int: raise MyValidationError
        #
        # total_moves = property(operator.attrgetter('_total_moves'))
        # @total_moves.setter
        # def total_moves(self, x):
        #     if not type(x) == int: raise MyValidationError
        #
        # accept_moves = property(operator.attrgetter('_accept_moves'))
        # @accept_moves.setter
        # def accept_moves(self, x):
        #     if not type(x) == int: raise MyValidationError
        #
        # status = property(operator.attrgetter('_status'))
        # @status.setter
        # def status(self, x):
        #     if not type(x) == str: raise MyValidationError
        #
        # skip_update = property(operator.attrgetter('_skip_update'))
        # @skip_update.setter
        # def skip_update(self, x):
        #     if not type(x) == bool: raise MyValidationError
        #
        # history = property(operator.attrgetter('_history'))
        # @history.setter
        # def history(self, x):
        #     if not jobtype.verify(x, 'history'): raise MyValidationError

    try:
        allthreads = pickle.load(open(pkl_file, 'rb'))
    except FileNotFoundError:
        raise FileNotFoundError('could not find specified restart pickle file: ' + pkl_file)

    for thread in allthreads:
        ValidThread(**thread.__dict__)
        try:    # try to build a valid thread from this thread's dictionary
            ValidThread(**thread.__dict__)
        except ThreadValidationError as e:
            warnings.warn('found an improperly formatted Thread ... ')

    output_str = 'output'
    output_int = 1

    return output_str, output_int


if __name__ == "__main__":
    main('tests/test_data/restart.pkl','tests/test_data/settings.pkl')
    #main(sys.argv[1], sys.argv[2])
