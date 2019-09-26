"""
Unit and regression test for utilities.py.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import sys
import pytest
import pytraj
import os
import glob
import filecmp
from atesa_v2 import utilities
from atesa_v2.configure import configure
from atesa_v2 import process

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa_v2/tests/test_temp'):
                os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
        except FileNotFoundError:
            pass
    
    def test_check_commit_fwd(self):
        """Tests check_commit using a dummy coordinate file and commitments defined to get result 'fwd'"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert utilities.check_commit('../test_data/test.rst7', settings) == 'fwd'

    def test_check_commit_bwd(self):
        """Tests check_commit using a dummy coordinate file and commitments defined to get result 'bwd'"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
        settings.commit_fwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert utilities.check_commit('../test_data/test.rst7', settings) == 'bwd'

    def test_check_commit_none(self):
        """Tests check_commit using a dummy coordinate file and commitments defined to get result ''"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert utilities.check_commit('../test_data/test.rst7', settings) == ''

    def test_check_commit_value_errors(self):
        """Tests check_commit using a dummy coordinate file and commitments defined to get ValueErrors"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        with pytest.raises(ValueError):
            utilities.check_commit('../test_data/test.rst7', settings)
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
        with pytest.raises(ValueError):
            utilities.check_commit('../test_data/test.rst7', settings)

    def test_get_cvs_no_qdot(self):
        """Tests get_cvs with a dummy coordinate file include_qdot = False"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = False
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        result = utilities.get_cvs('../test_data/test.rst7', settings)
        assert len(result.split(' ')) == 3   # includes empty '' at end
        assert float(result.split(' ')[0]) == pytest.approx(1.09, 1e-3)
        assert float(result.split(' ')[1]) == pytest.approx(1.69, 1e-3)

    def test_get_cvs_with_qdot_broken(self):
        """Tests get_cvs with a dummy coordinate file include_qdot = True and a coordinate file without velocities"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = True
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        with pytest.raises(IndexError):
            result = utilities.get_cvs('../test_data/test.rst7', settings)

    def test_get_cvs_with_qdot(self):
        """Tests get_cvs with a dummy coordinate file include_qdot = True and a coordinate file with velocities"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = True
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        result = utilities.get_cvs('../test_data/test_velocities.rst7', settings)
        assert len(result.split(' ')) == 5   # includes empty '' at end
        assert float(result.split(' ')[0]) == pytest.approx(1.089, 1e-3)
        assert float(result.split(' ')[1]) == pytest.approx(1.728, 1e-3)
        assert float(result.split(' ')[2]) == pytest.approx(-0.252, 1e-2)
        assert float(result.split(' ')[3]) == pytest.approx(0.282, 1e-2)

    def test_rev_vels(self):
        """Tests rev_vels by comparing to a known-correct file"""
        rev_file = utilities.rev_vels('../test_data/test_velocities.rst7')
        assert filecmp.cmp('../test_data/test_velocities_reversed.rst7', rev_file)
        os.remove(rev_file)

    def test_evaluate_rc(self):
        """Tests evaluate_rc"""
        cv_list = [1, 0.3, -3, 1.5002]  # various data types and precisions
        rc_defn = '2*CV0 + 0.3*CV1 + CV2/0.4 - numpy.sin(CV3)'
        assert utilities.evaluate_rc(rc_defn, cv_list) == pytest.approx(-6.408, 1E-2)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
