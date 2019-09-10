"""
Unit and regression test for process.py.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import pytest
import sys
from ..configure import configure
from ..process import process
import pytraj
import os
import glob

class Tests(object):
    def test_process_already_terminated(self):
        """Tests process.py for a thread with terminated = True"""
        settings = configure('atesa_v2/data/atesa.config')
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].terminated = True
        running = [allthreads[0]]
        assert process(allthreads[0], running, settings) == []

    def test_process_to_be_terminated(self):
        """Tests process.py for a thread with terminated = False and current_type = 'terminate'"""
        settings = configure('atesa_v2/data/atesa.config')
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].terminated = False
        allthreads[0].current_type = 'terminate'
        running = [allthreads[0]]
        assert process(allthreads[0], running, settings) == []

    def test_process_to_submit(self):
        """Tests process.py for a thread that should be submitted"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.DEBUG = True
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].terminated = False
        allthreads[0].current_type = []
        running = []
        try:
            os.chdir('atesa_v2/tests/test_temp')
        except FileNotFoundError:
            os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
        assert process(allthreads[0], running, settings) == [allthreads[0]]
        assert allthreads[0].jobids == ['123456']
        assert allthreads[0].current_type == ['init']
        assert allthreads[0].current_name == ['init']
        os.chdir('../../..')    # to go back to the original testing directory

    @classmethod
    def teardown_class(cls):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
