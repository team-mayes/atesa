"""
Unit and regression test for utilities.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import pytest
import pytraj
import os
import glob
import filecmp
import shutil
import pickle
from atesa import utilities
from atesa.configure import configure
from atesa import process
from atesa import main

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
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

    def test_check_commit_traj(self):
        """Tests check_commit using a dummy trajectory file and commitments defined to get result 'fwd'"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert utilities.check_commit('../test_data/test.nc', settings) == 'fwd'

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
        """Tests get_cvs with a dummy coordinate file and include_qdot = False"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = False
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        result = utilities.get_cvs('../test_data/test.rst7', settings)
        assert len(result.split(' ')) == 2
        assert float(result.split(' ')[0]) == pytest.approx(1.09, 1e-3)
        assert float(result.split(' ')[1]) == pytest.approx(1.69, 1e-3)

    def test_get_cvs_with_qdot_broken(self):
        """Tests get_cvs with a dummy coordinate file, include_qdot = True, and a coordinate file without velocities"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = True
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        with pytest.raises(IndexError):
            result = utilities.get_cvs('../test_data/test.rst7', settings)

    def test_get_cvs_with_qdot(self):
        """Tests get_cvs with a dummy coordinate file, include_qdot = True, and a coordinate file with velocities"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = True
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        result = utilities.get_cvs('../test_data/test_velocities.rst7', settings)
        assert len(result.split(' ')) == 4
        assert float(result.split(' ')[0]) == pytest.approx(1.089, 1e-3)
        assert float(result.split(' ')[1]) == pytest.approx(1.728, 1e-3)
        assert float(result.split(' ')[2]) == pytest.approx(-0.252, 1e-2)
        assert float(result.split(' ')[3]) == pytest.approx(0.282, 1e-2)

    def test_get_cvs_no_qdot_reduce(self):
        """Tests get_cvs with a dummy coordinate file, include_qdot = False, and reduce = True"""
        settings = configure('../../data/atesa.config')
        settings.include_qdot = False
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @3\')[0]', 'pytraj.distance(traj, \'@2 @4\')[0]']
        with pytest.raises(FileNotFoundError):  # no as.out present yet
            result = utilities.get_cvs('../test_data/test.rst7', settings, reduce=True)
        shutil.copy('../test_data/as.out', 'as.out')    # copy in the necessary file
        result = utilities.get_cvs('../test_data/test.rst7', settings, reduce=True)
        assert len(result.split(' ')) == 2
        assert float(result.split(' ')[0]) == pytest.approx(0.571, 1e-2)
        assert float(result.split(' ')[1]) == pytest.approx(0.326, 1e-2)

    def test_rev_vels(self):
        """Tests rev_vels by comparing to a known-correct file"""
        rev_file = utilities.rev_vels('../test_data/test_velocities.rst7')
        assert filecmp.cmp('../test_data/test_velocities_reversed.rst7', rev_file)
        os.remove(rev_file)

    def test_evaluate_rc(self):
        """Tests evaluate_rc"""
        cv_list = [1, 0.3, -3, 1.5002]  # various data types and precisions
        rc_defn = '2*CV1 + 0.3*CV2 + CV3/0.4 - numpy.sin(CV4)'
        assert utilities.evaluate_rc(rc_defn, cv_list) == pytest.approx(-6.408, 1E-2)

    def test_resample(self):
        """Tests resample"""
        settings = configure('../../data/atesa.config')
        with pytest.raises(FileNotFoundError):      # run without making restart.pkl
            utilities.resample(settings)
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test.rst7', '../test_data/test_two_init.rst7']
        settings.rc_definition = '1.00 + 2.34*CV1 - 5.67*CV2'
        settings.rc_reduced_cvs = False
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.include_qdot = False
        settings.topology = '../test_data/test.prmtop'
        allthreads = main.init_threads(settings)
        allthreads[0].history.init_coords = [['../test_data/test_velocities_init.rst7', '../test_data/test_velocities_init_bwd.rst7'],
                                             ['../test_data/test_two_init.rst7', '../test_data/test_two_init_bwd.rst7']]
        allthreads[0].history.prod_results = [['fwd', 'bwd'], ['bwd', 'bwd']]
        allthreads[0].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc'], ['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[1].history.init_coords = [['../test_data/test_velocities_init.rst7', '../test_data/test_velocities_init_bwd.rst7'],
                                             ['../test_data/test_two_init.rst7', '../test_data/test_two_init_bwd.rst7']]
        allthreads[1].history.prod_results = [['fwd', 'fwd'], ['bwd', 'fwd']]
        allthreads[1].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc'], ['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[0].history.timestamps = [1, 3]
        allthreads[1].history.timestamps = [2, 4]
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)  # file now exists
        utilities.resample(settings, full_cvs=True)
        assert os.path.exists('as_raw.out')
        # assert os.path.exists('as_decorr.out')    # todo: reinstitute this test once I have a good restart.pkl file to test with; todo: figure out what I meant by this
        raw_lines = open('as_raw.out', 'r').readlines()
        assert len(raw_lines) == 4
        assert not False in [raw_lines[0][0] == 'B', raw_lines[1][0] == 'A', raw_lines[2][0] == 'B', raw_lines[3][0] == 'A']
        new_lines = open('as_full_cvs.out', 'r').readlines()
        cmp_lines = open('../test_data/as_full_cvs.out', 'r').readlines()
        for line_index in range(len(new_lines)):
            for item_index in range(len(new_lines[line_index].split())):
                assert pytest.approx(float(new_lines[line_index].split()[item_index]), float(cmp_lines[line_index].split()[item_index]), 1e-3)

    def test_interpret_cv(self):
        """Test interpret_cv for a distance, angle, dihedral, and difference of distances for both pytraj and mdtraj"""
        settings = configure('../../data/atesa.config')
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]',
                        'pytraj.angle(traj, \'@1 @2 @3\')[0]',
                        'pytraj.dihedral(traj, \'@1 @2 @3 @4\')[0]',
                        'pytraj.distance(traj, \'@1 @2\')[0] - pytraj.distance(traj, \'@3 @4\')[0]',
                        'mdtraj.compute_distances(mtraj, numpy.array([[83, 84]]))[0][0] * 10',
                        'mdtraj.compute_angles(mtraj, numpy.array([[83, 85, 86]]))[0][0] * 180 / numpy.pi',
                        'mdtraj.compute_dihedrals(mtraj, numpy.array([[84, 83, 85, 86]]))[0][0] * 180 / numpy.pi',
                        '(mdtraj.compute_distances(mtraj, numpy.array([[101, 103]]))[0][0] * 10) - (mdtraj.compute_distances(mtraj, numpy.array([[101, 104]]))[0][0] * 10)']
        assert utilities.interpret_cv(1, settings) == (['1', '2'], 'distance', 2)
        assert utilities.interpret_cv(2, settings) == (['1', '2', '3'], 'angle', 3)
        assert utilities.interpret_cv(3, settings) == (['1', '2', '3', '4'], 'dihedral', 4)
        assert utilities.interpret_cv(4, settings) == (['1', '2', '3', '4'], 'diffdistance', 4)
        assert utilities.interpret_cv(5, settings) == (['84', '85'], 'distance', 2)
        assert utilities.interpret_cv(6, settings) == (['84', '86', '87'], 'angle', 3)
        assert utilities.interpret_cv(7, settings) == (['85', '84', '86', '87'], 'dihedral', 4)
        assert utilities.interpret_cv(8, settings) == (['102', '104', '102', '105'], 'diffdistance', 4)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
