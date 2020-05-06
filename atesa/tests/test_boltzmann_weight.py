"""
Unit and regression test for boltzmann_weight.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import os
import shutil
import pytest
import pickle
import subprocess
import glob
from atesa import boltzmann_weight
from atesa.configure import configure

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_main(self):
        """Tests boltzmann_weight using example EPS output in test_data"""
        kwargs = {'i': '../test_data/eps_example.out', 'o': 'fep_test.out', 't': 350, 'n': 4, 'b': -1, 'c': 1, 'noplot': True, 'e': 0.8, 'slope_only': False}
        # There's a necessary random component here that can potentially cause a RuntimError to be raised, but the odds
        # of it failing three times in a row without there actually being something wrong are extremely small, so this
        # is the best I can think to do...
        try:
            boltzmann_weight.main(**kwargs)
        except RuntimeError:
            try:
                boltzmann_weight.main(**kwargs)
            except RuntimeError:
                boltzmann_weight.main(**kwargs)
        compare_lines = open('../test_data/fep_example.out', 'r').readlines()
        line_index = 0
        for line in open('fep_test.out', 'r').readlines():
            # Compare first and second entries in each line (third is bootstrapped and so will vary)
            assert line.split()[0] == compare_lines[line_index].split()[0]
            assert line.split()[1] == compare_lines[line_index].split()[1]
            line_index += 1

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
