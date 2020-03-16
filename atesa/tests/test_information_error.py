"""
Unit and regression test for information_error.py.
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
from atesa import information_error
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
        """Tests information_error.main using sham shooting points in test_data"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.include_qdot = False
        settings.DEBUG = True   # to skip resampling and substitute our own as_decorr.out file
        shutil.copy('../test_data/test_velocities_init.rst7', '../test_temp')
        shutil.copy('../test_data/test_two_init.rst7', '../test_temp')
        shutil.copy('../test_data/as.out', '../test_temp')
        shutil.copy('../test_data/as_big.out', 'as_decorr_6054.out')

        # First test without extant settings.pkl file
        with pytest.raises(FileNotFoundError):
            information_error.main()

        settings.__dict__.pop('env')    # env attribute is not picklable
        pickle.dump(settings, open('settings.pkl', 'wb'), protocol=2)   # main will look for this file to load in settings

        shutil.copy('../test_data/as_big.out', 'as_decorr_6054.out')

        information_error.main()
        assert os.path.exists('info_err.out')
        assert float(open('info_err.out', 'r').readlines()[0].split(' ')[0]) == 6054
        assert float('%.2f' % float(open('info_err.out', 'r').readlines()[0].split(' ')[1])) == pytest.approx(0.16, 1E-2)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
