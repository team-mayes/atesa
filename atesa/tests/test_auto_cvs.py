"""
Unit and regression test for auto_cvs.py.
"""

# Import package, test suite, and other packages as needed
import atesa
import sys
import os
import shutil
import pytest
import glob
import filecmp
import argparse
import mdtraj
import pytraj
import numpy
import scipy
from atesa import auto_cvs
from atesa.configure import configure

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_main_mdtraj(self):
        """Tests auto_cvs with mdtraj"""
        settings = argparse.Namespace()
        settings.auto_cvs_type = 'mdtraj'
        settings.topology = '../test_data/test.prmtop'
        settings.initial_coordinates = ['../test_data/test.rst7']
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.commit_fwd = ([101, 102, 102], [103, 104, 105], [1.5, 2.0, 1.0], ['lt', 'gt', 'lt'])
        settings.commit_bwd = ([101, 102, 102], [103, 104, 105], [2.0, 1.5, 1.8], ['gt', 'lt', 'gt'])
        settings.working_directory = './'
        settings.auto_cvs_radius = 5
        settings.auto_cvs_exclude_water = True

        mtraj = mdtraj.load(settings.initial_coordinates[0], top=settings.topology)  # for testing third assert

        result = auto_cvs.main(settings)
        compare = open('../test_data/cvs_mdtraj.txt', 'r').readlines()[1:]

        assert len(result) == len(compare)
        assert all([item in [line.split('; ')[1].strip('\n') for line in compare] for item in result])
        assert pytest.approx(eval(compare[0].split('; ')[1]), 1.09, 1E-2)

    def test_main_pytraj(self):
        """Tests auto_cvs with pytraj"""
        settings = argparse.Namespace()
        settings.auto_cvs_type = 'pytraj'
        settings.topology = '../test_data/test.prmtop'
        settings.initial_coordinates = ['../test_data/test.rst7']
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.commit_fwd = ([101, 102, 102], [103, 104, 105], [1.5, 2.0, 1.0], ['lt', 'gt', 'lt'])
        settings.commit_bwd = ([101, 102, 102], [103, 104, 105], [2.0, 1.5, 1.8], ['gt', 'lt', 'gt'])
        settings.working_directory = './'
        settings.auto_cvs_radius = 5
        settings.auto_cvs_exclude_water = True

        traj = pytraj.load(settings.initial_coordinates[0], settings.topology)  # for testing third assert

        result = auto_cvs.main(settings)
        compare = open('../test_data/cvs_pytraj.txt', 'r').readlines()[1:]

        assert len(result) == len(compare)
        assert all([item in [line.split('; ')[1].strip('\n') for line in compare] for item in result])
        assert pytest.approx(eval(compare[0].split('; ')[1]), 1.09, 1E-2)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)

