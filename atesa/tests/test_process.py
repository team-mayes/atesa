"""
Unit and regression test for process.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import pytest
import sys
import pytraj
import os
import glob
from atesa.configure import configure
from atesa.process import process
from atesa import main

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_process_to_be_terminated(self):
        """Tests process.py for a thread with terminated = True"""
        settings = configure('../../data/atesa.config')
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].terminated = True
        running = [allthreads[0]]
        process(allthreads[0], running, settings)
        assert running == []

    def test_process_to_submit(self):
        """Tests process.py for a thread that should be submitted"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.DEBUG = True
        allthreads = main.init_threads(settings)
        allthreads[0].terminated = False
        allthreads[0].current_type = []
        running = []
        assert process(allthreads[0], running, settings) == [allthreads[0]]
        assert allthreads[0].jobids == ['123456']
        assert allthreads[0].current_type == ['init']
        assert allthreads[0].current_name == ['init']

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
