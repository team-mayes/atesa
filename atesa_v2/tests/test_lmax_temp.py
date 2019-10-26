"""
Unit and regression test for lmax.py.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import sys
import pytest
import os
import shutil
import glob
import argparse
from atesa_v2 import lmax_temp

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa_v2/tests/test_temp'):
                os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
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
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
