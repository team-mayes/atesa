"""
Unit and regression test for utilities.py.
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
from atesa import rc_eval
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
        """Tests rc_eval.main using sham shooting points in test_data"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.include_qdot = False
        shutil.copy('../test_data/test_velocities_init.rst7', '../test_temp')
        shutil.copy('../test_data/test_two_init.rst7', '../test_temp')
        shutil.copy('../test_data/as.out', '../test_temp')

        # First test without extant settings.pkl file
        with pytest.raises(FileNotFoundError):
            rc_eval.main('../test_temp/', '3*CV1 - 0.3*CV2', '../test_data/as.out')
        process = subprocess.Popen('rc_eval.py ../test_temp/ 3*CV0-0.3*CV1 ../test_data/as.out', stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True, shell=True)
        output = process.stdout.read().decode()
        assert 'FileNotFoundError' in output    # this section just to ensure it works when called both ways

        settings.__dict__.pop('env')    # env attribute is not picklable
        pickle.dump(settings, open('settings.pkl', 'wb'), protocol=2)   # main will look for this file to load in settings
        # shutil.move('settings.pkl', '../test_data/settings.pkl')
        rc_eval.main('../test_temp/', '3*CV1 - 0.3*CV2', '../test_data/as.out')
        assert os.path.exists('rc.out')
        lines = open('rc.out', 'r').readlines()
        for i in range(len(lines)):     # assert that the rc.out file is properly sorted
            try:
                assert abs(float(lines[i].split(': ')[1])) < abs(float(lines[i+1].split(': ')[1]))
            except IndexError:
                pass

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
