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
from atesa import lmax

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
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
            lmax.two_line_test_func(results, False)
        # Next, data designed not to pass due to too mild a slope change
        results = data_converter([-100, -110, -120, -129, -138])
        assert lmax.two_line_test_func(results, False) == -1
        # Data designed to pass at index 2 because it's just garbage (very much a fringe test)
        results = data_converter([-100, -110, -120, -122, -124])
        assert lmax.two_line_test_func(results, False) == 2

    def test_objective_function(self):
        """Tests objective_function using dummy data"""
        # Assemble data with known outcome
        A_list = [[1, 2], [1.1, 2.1]]
        B_list = [[0.5, 2.5], [0.3, 2.6]]
        params = [2.25, 1.12, -0.94]
        assert lmax.objective_function(params, A_list, B_list) == pytest.approx(8.99, 1E-2)

    def test_main_non_existent_input(self):
        """Tests main using an input file name that does not exist"""
        kwargs = {'two_line_test': True, 'i': ['definitely_not_a_real_file_who_would_name_a_file_this.txt'], 'k': [0], 'f': [], 'q': ['absent'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        with pytest.raises(FileNotFoundError):
            lmax.main(**kwargs)

    def test_main_improper_input(self):
        """Tests main using an input file that does exist, but is improperly formatted"""
        kwargs = {'two_line_test': True, 'i': ['../test_data/test.rst7'], 'k': [0], 'f': [], 'q': ['absent'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        with pytest.raises(RuntimeError):
            lmax.main(**kwargs)

    def test_main_two_line_test(self):
        """Tests main using two_line_test"""
        kwargs = {'two_line_test': True, 'i': ['../test_data/as.out'], 'k': [0], 'f': [], 'q': ['absent'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        lmax.main(**kwargs)
        results = open('lmax.out', 'r').readlines()[1].split(' + ')
        assert float(results[0].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '')) == pytest.approx(5.47, 1E-2)
        assert float(results[1].replace('*CV1', '')) == pytest.approx(-22.116, 1E-2)
        assert float(results[2].replace('*CV4', '')) == pytest.approx(5.483, 1E-2)
        assert len(results) == 3

    def test_main_running(self):
        """Tests main using running"""
        kwargs = {'two_line_test': False, 'i': ['../test_data/as.out'], 'k': [0], 'f': [], 'q': ['absent'], 'r': [2], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        lmax.main(**kwargs)
        results = open('lmax.out', 'r').readlines()[1].split(' + ')
        assert float(results[0].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '')) == pytest.approx(5.47, 1E-2)
        assert float(results[1].replace('*CV1', '')) == pytest.approx(-22.116, 1E-2)
        assert float(results[2].replace('*CV4', '')) == pytest.approx(5.483, 1E-2)
        assert len(results) == 3

    # This test deprecated for now because it's unimportant and causing CI issues due to platform differences in which
    # CVs get included (since the log likelihood is so tiny in each case for this tiny as.out file)
    # def test_main_k(self):
    #     """Tests main using k and not f"""
    #     kwargs = {'two_line_test': False, 'i': ['../test_data/as.out'], 'k': [2], 'f': [], 'q': ['absent'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55]}
    #     lmax.main(**kwargs)
    #     results = open('lmax.out', 'r').readlines()[1].split(' + ')
    #     assert float(results[0].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '')) == pytest.approx(-6.494, 1E-2)
    #     assert float(results[1].replace('*CV2', '')) == pytest.approx(7.959, 1E-2)
    #     assert float(results[2].replace('*CV5', '')) == pytest.approx(2.247, 1E-2)
    #     assert len(results) == 3

    def test_main_k_and_f(self):
        """Tests main using k and f"""
        kwargs = {'two_line_test': False, 'i': ['../test_data/as.out'], 'k': [2], 'f': [1, 4], 'q': ['absent'], 'r': [0], 'o': ['lmax.out'], 'quiet': False, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        lmax.main(**kwargs)
        results = open('lmax.out', 'r').readlines()[1].split(' + ')
        assert float(results[0].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '')) == pytest.approx(5.47, 1E-2)
        assert float(results[1].replace('*CV1', '')) == pytest.approx(-22.116, 1E-2)
        assert float(results[2].replace('*CV4', '')) == pytest.approx(5.483, 1E-2)
        assert len(results) == 3

    def test_main_k_f_and_q(self):
        """Tests main using k, f, and q"""
        kwargs = {'two_line_test': False, 'i': ['../test_data/as.out'], 'k': [2], 'f': [1, 4], 'q': ['present'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        with pytest.raises(RuntimeError):   # this as.out file has an uneven number of CVs
            lmax.main(**kwargs)
        kwargs = {'two_line_test': False, 'i': ['as.out'], 'k': [2], 'f': [1, 2], 'q': ['present'], 'r': [0], 'o': ['lmax.out'], 'quiet': True, 'plots': False, 'two_line_threshold': [0.55], 's': [], 'hist_bins': [10]}
        shutil.copy('../test_data/as.out', 'as_temp.out')   # need to make a new as.out file with an even number of CVs
        open('as.out', 'w').close()
        for line in open('as_temp.out', 'r').readlines():
            line = line.strip('\n') + ' 0.01\n'
            open('as.out', 'a').write(line)
        lmax.main(**kwargs)
        results = open('lmax.out', 'r').readlines()[1].split(' + ')
        assert float(results[0].replace('The optimized reaction coordinate (with CVs indexed from 1) is: ', '')) == pytest.approx(-0.654, 1E-2)
        assert float(results[1].replace('*CV1', '')) == pytest.approx(-6.279, 1E-2)
        assert float(results[2].replace('*CV2', '')) == pytest.approx(11.096, 1E-2)
        assert len(results) == 3

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
