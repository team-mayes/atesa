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
from atesa_v2 import lmax

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa_v2/tests/test_temp'):
                os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_two_line_test(self):
        """Tests two_line_test using dummy data"""
        def data_converter(my_list):
            # Helper function to convert data types quickly
            out = []
            for item in my_list:
                temp = argparse.Namespace()
                temp.fun = item
                out.append(temp)

            return out

        # First assemble data designed not to pass due to insufficient data:
        results = data_converter([-100, -110, -120, -130])
        with pytest.raises(RuntimeError):
            lmax.two_line_test(results)
        # Next, data designed not to pass due to too mild a slope change
        results = data_converter([-100, -110, -120, -125, -130])
        assert lmax.two_line_test(results) == -1
        # Data designed to pass at index 2 because it's just garbage (very much a fringe test)
        results = data_converter([-100, -110, -120, -122, -124])
        assert lmax.two_line_test(results) == 2

    def test_objective_function(self):
        """Tests objective_function using dummy data"""
        # Assemble data with known outcome
        A_list = [[1, 2], [1.1, 2.1]]
        B_list = [[0.5, 2.5], [0.3, 2.6]]
        params = [2.25, 1.12, -0.94]
        assert lmax.objective_function(params, A_list, B_list) == pytest.approx(8.99, 1E-2)

    def test_main_non_existent_input(self):
        """Tests main using an input file name that does not exist"""
        kwargs = {'automagic': True, 'i': ['definitely_not_a_real_file_who_would_name_a_file_this.txt']}
        with pytest.raises(FileNotFoundError):
            lmax.main(**kwargs)

    def test_main_improper_input(self):
        """Tests main using an input file that does exist, but is improperly formatted"""
        kwargs = {'automagic': True, 'i': ['../test_data/test.rst7']}
        with pytest.raises(RuntimeError):
            lmax.main(**kwargs)

    def test_main_automagic(self):
        """Tests main using automagic"""
        kwargs = {'automagic': True, 'i': ['../test_data/as.out'], 'k': [0], 'f': [], 'q': [False], 'r': [0], 'o': ['lmax.out'], 'quiet': True}
        lmax.main(**kwargs)
        assert '5.470 + -22.116*CV1 + 5.483*CV4' in open('lmax.out', 'r').readlines()[1]

    def test_main_k(self):
        """Tests main using k and not f"""
        kwargs = {'automagic': False, 'i': ['../test_data/as.out'], 'k': [2], 'f': [], 'q': [False], 'r': [0], 'o': ['lmax.out'], 'quiet': True}
        lmax.main(**kwargs)
        assert '-6.494 + 7.959*CV2 + 2.247*CV5' in open('lmax.out', 'r').readlines()[1]     # different from automagic result because this doesn't use running

    def test_main_k_and_f(self):
        """Tests main using k and f"""
        kwargs = {'automagic': False, 'i': ['../test_data/as.out'], 'k': [2], 'f': [1, 4], 'q': [False], 'r': [0], 'o': ['lmax.out'], 'quiet': True}
        lmax.main(**kwargs)
        assert '5.470 + -22.116*CV1 + 5.483*CV4' in open('lmax.out', 'r').readlines()[1]

    def test_main_k_f_and_q(self):
        """Tests main using k, f, and q"""
        kwargs = {'automagic': False, 'i': ['as.out'], 'k': [2], 'f': [1], 'q': [True], 'r': [0], 'o': ['lmax.out'], 'quiet': True}
        shutil.copy('../test_data/as.out', 'as_temp.out')   # need to make a new as.out file with an even number of CVs
        open('as.out', 'w').close()
        for line in open('as_temp.out', 'r').readlines():
            line = line.strip('\n') + ' 0.01\n'
            open('as.out', 'a').write(line)
        lmax.main(**kwargs)
        assert '-0.654 + -6.279*CV1 + 11.096*CV2' in open('lmax.out', 'r').readlines()[1]

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
