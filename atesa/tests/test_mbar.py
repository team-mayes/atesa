"""
Unit and regression test for mbar.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import pytest
import os
import shutil
import glob
import argparse
import filecmp
from atesa import mbar

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_basic_settings(self):
        """Tests mbar.py with basic (default) settings, except quiet = True"""
        for file in glob.glob('../test_data/rcwin_*_us.dat'):
            shutil.copy(file, './')
        kwargs = {'t': [300], 'k': [50], 'o': ['mbar.out'], 'min_data': [0], 'ignore': [1], 'decorr': False, 'rc_min': [''], 'rc_max': [''], 'quiet': True, 'i': ['./']}
        mbar.main(**kwargs)
        assert filecmp.cmp('../test_data/mbar.out', 'mbar.out')

    def test_no_data(self):
        """Tests mbar.py with no valid files in the working directory"""
        kwargs = {'t': [300], 'k': [50], 'o': ['mbar.out'], 'min_data': [0], 'ignore': [1], 'decorr': False,
                  'rc_min': [''], 'rc_max': [''], 'quiet': True, 'i': ['./']}
        with pytest.raises(RuntimeError):
            mbar.main(**kwargs)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
