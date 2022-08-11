"""
Unit and regression test for resampling.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import os
import shutil
import pytest
import pickle
import glob
import filecmp
import argparse
from atesa.configure import configure
from atesa import main
from atesa.resampling import resample_committor_analysis
from atesa.resampling import resample_umbrella_sampling

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_resample_committor_analysis(self):
        """Tests resample_committor_analysis"""
        # First spoof settings file and trajectories in working directory
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.resample = False   # explicitly false in spoofed settings file
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        settings.committor_analysis_use_rc_out = False
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.5, 2.0], ['lt', 'gt']]    # doesn't match test.nc
        settings.commit_bwd = [[1, 2], [3, 4], [2.0, 1.5], ['gt', 'lt']]    # doesn't match test.nc
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod', 'prod']
        shutil.copy('../test_data/test.nc', 'test.rst7_0_0.nc')
        shutil.copy('../test_data/test.nc', 'test.rst7_0_1.nc')
        shutil.copy('../test_data/test.nc', 'test.rst7_0_2.nc')
        shutil.copy('../test_data/test.rst7', 'test.rst7')
        allthreads[0].history.prod_trajs = ['test.nc', 'test.nc', 'test.nc']
        pickle.dump(settings, open('settings.pkl', 'wb'))

        # Now make settings for resampling
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.resample = True
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        settings.committor_analysis_use_rc_out = False
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.1, 1.6], ['lt', 'gt']]    # DOES match test.nc
        settings.commit_bwd = [[1, 2], [3, 4], [1.6, 1.1], ['gt', 'lt']]    # doesn't match test.nc

        # Finally, call resample
        resample_committor_analysis(settings)

        allthreads = pickle.load(open('restart.pkl', 'rb')) # reload resampled allthreads
        assert os.path.exists('committor_analysis.out')
        assert len(allthreads) == 1     # just one thread
        assert allthreads[0].history.prod_results == ['fwd', 'fwd', 'fwd']

    def test_resample_umbrella_sampling(self):
        """Tests resample_umbrella_sampling"""
        # First spoof trajectories in working directory
        shutil.copy('../test_data/test.nc', 'test.rst7_0_0.nc')
        shutil.copy('../test_data/test.nc', 'test.rst7_0_1.nc')
        shutil.copy('../test_data/test.nc', 'test.rst7_0_2.nc')

        # Now make settings for resampling
        settings = configure('../../data/atesa.config')
        settings.job_type = 'umbrella_sampling'
        settings.resample = True
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.rc_definition = '1.00 + 2.34*CV1 - 0.67*CV2'

        # Finally, call resample
        resample_umbrella_sampling(settings)

        assert os.path.exists('us_full_cvs.out')

        lines = open('us_full_cvs.out', 'r').readlines()

        assert len(lines) == 6
        assert len(lines[0].split()) == 2

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)

