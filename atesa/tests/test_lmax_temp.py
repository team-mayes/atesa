"""
Unit and regression test for lmax.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import pytest
import os
import shutil
import glob
import argparse
from atesa import lmax_temp

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    # def test_main(self):
    #     """Tests main"""
    #     err = lmax_temp.main()
    #     print(err)
    #     assert err == pytest.approx(0.209, 1E-2)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
