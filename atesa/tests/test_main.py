"""
Unit and regression test for the atesa_v2 package.
"""

# Import package, test suite, and other packages as needed
import atesa
import pytest
import sys
import pytraj
import os
import glob
import shutil
import pickle
from atesa.configure import configure
from atesa import main

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass
    
    def test_atesa_v2_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "atesa" in sys.modules

    def test_main(self):
        """Tests main.main with DEBUG = True to skip actual job submission/monitoring"""
        settings = configure('../../data/atesa.config')
        settings.initial_coordinates = ['../test_data/test.rst7']
        settings.topology = '../test_data/test.prmtop'
        settings.degeneracy = 2
        settings.overwrite = True
        settings.resample = False
        settings.DEBUG = True
        assert main.main(settings) == 'ATESA run exiting normally'

    def test_configure_broken(self):
        """Tests configure.py with a non-existent file"""
        with pytest.raises(FileNotFoundError):
            settings = configure('../../data/ates.config')

    def test_configure_directory(self):
        """Tests configure.py behavior in correcting an improperly formatted directory"""
        shutil.copy('../../data/atesa.config', 'atesa.config')  # todo: make a dedicated 'test' copy of atesa.config for tests to use.
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
        assert allthreads[0].topology == 'test.prmtop'
        assert allthreads[1].topology == 'test.prmtop'

    def test_init_threads_restart(self):
        """Tests successful initialization of restarted threads"""
        settings = configure('../../data/atesa.config')
        # First, make restart.pkl
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test.rst7', '../test_data/test_two_init.rst7']
        settings.rc_definition = '1.00 + 2.34*CV0 - 5.67*CV1'
        settings.rc_reduced_cvs = False
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.include_qdot = False
        settings.topology = '../test_data/test.prmtop'
        allthreads = main.init_threads(settings)
        allthreads[0].history.init_coords = [['../test_data/test_velocities_init.rst7', '../test_data/test_velocities_init_bwd.rst7'],
                                             ['../test_data/test_two_init.rst7', '../test_data/test_two_init_bwd.rst7']]
        allthreads[0].history.prod_results = [['fwd', 'bwd'], ['bwd', 'bwd']]
        allthreads[1].history.init_coords = [['../test_data/test_velocities_init.rst7', '../test_data/test_velocities_init_bwd.rst7'],
                                             ['../test_data/test_two_init.rst7', '../test_data/test_two_init_bwd.rst7']]
        allthreads[1].history.prod_results = [['bwd', 'bwd'], ['fwd', 'bwd']]
        allthreads[0].history.timestamps = [1, 3]
        allthreads[1].history.timestamps = [2, 4]
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)
        # Then, reset settings, set restart = True, and try again
        settings = configure('../../data/atesa.config')
        settings.restart = True
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        assert len(allthreads) == 2
        assert allthreads[0].history.prod_results == [['fwd', 'bwd'], ['bwd', 'bwd']]
        assert allthreads[1].history.prod_results == [['bwd', 'bwd'], ['fwd', 'bwd']]
        assert allthreads[0].topology == 'test.prmtop'
        assert allthreads[1].topology == 'test.prmtop'

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

    def test_import_cvs(self):
        """Tests importing CVs and commitment definitions from a given settings.pkl file"""
        #shutil.copy('../../data/atesa.config', 'temp.config')
        open('temp.config', 'w').write('job_type = \'committor_analysis\'')
        open('temp.config', 'a').write('\nbatch_system = \'slurm\'')
        open('temp.config', 'a').write('\nrestart = False')
        open('temp.config', 'a').write('\ntopology = \'../test_data/test.prmtop\'')
        open('temp.config', 'a').write('\nworking_directory = \'./\'')
        open('temp.config', 'a').write('\noverwrite = False')
        open('temp.config', 'a').write('\nas_settings_file = \'../test_data/settings.pkl\'')
        settings = configure('temp.config')
        assert settings.cvs == ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        assert settings.commit_fwd == [[1, 2], [3, 4], [1.5, 2.0], ['lt', 'gt']]
        assert settings.commit_bwd == [[1, 2], [3, 4], [2.0, 1.5], ['gt', 'lt']]

    @classmethod
    def teardown_method(self, method):
        "Runs at end of each method"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)
