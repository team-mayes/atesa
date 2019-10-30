"""
Unit and regression test for the atesa_v2 package.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import pytest
import sys
import pytraj
import os
import glob
import shutil
import pickle
from atesa_v2.configure import configure
from atesa_v2 import main

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa_v2/tests/test_temp'):
                os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
        except FileNotFoundError:
            pass
    
    def test_atesa_v2_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "atesa_v2" in sys.modules

    def test_main(self):
        """Tests main.main with no initial coordinates given"""
        settings = configure('../../data/atesa.config')
        settings.initial_coordinates = []
        settings.topology = '../test_data/test.prmtop'
        settings.overwrite = True
        settings.resample = False
        assert main.main(settings) == 'ATESA run exiting normally (all threads ended individually)'

    def test_configure_broken(self):
        """Tests configure.py with a non-existent file"""
        with pytest.raises(FileNotFoundError):
            settings = configure('../../data/ates.config')

    def test_configure_directory(self):
        """Tests configure.py behavior in correcting an improperly formatted directory"""
        shutil.copy('../../data/atesa.config', 'atesa.config')  # todo: before publication, make a dedicated 'test' copy of atesa.config for tests to use.
        config_lines = open('atesa.config', 'r').readlines()
        for line in config_lines:
            if 'working_directory' in line:
                line = 'working_directory = \'/foo/bar/\'\n'
            open('atesa_temp.config', 'a').write(line)
        settings = configure('atesa_temp.config')
        assert settings.working_directory == '/foo/bar'

    def test_init_threads_new(self):
        """Tests successful initialization of new threads"""
        settings = configure('../../data/atesa.config')
        settings.initial_coordinates = ['../test_data/test.rst7', '../test_data/test_two_init.rst7']
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        assert len(allthreads) == 2
        assert allthreads[0].history.init_inpcrd == ['test.rst7']
        assert allthreads[1].history.init_inpcrd == ['test_two_init.rst7']
        assert allthreads[0].topology == 'topology.prmtop'
        assert allthreads[1].topology == 'topology.prmtop'

    def test_thread_get_frame_amber(self):
        """Tests thread.get_frame method with md_engine = 'amber' and frame = -1"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.md_engine = 'amber'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].jobids.append('01234')
        allthreads[0].history.prod_trajs.append(['../test_data/test.nc'])
        compare_traj = pytraj.iterload('../test_data/test.rst7', allthreads[0].topology)
        test_frame = allthreads[0].get_frame(allthreads[0].history.prod_trajs[0][0], -1, settings)
        query_traj = pytraj.iterload(test_frame, allthreads[0].topology)
        assert query_traj.n_frames == compare_traj.n_frames
        assert pytraj.center_of_mass(query_traj) == pytest.approx(pytraj.center_of_mass(compare_traj), 1e-3)

    @classmethod
    def teardown_method(self, method):
        "Runs at end of each method"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
