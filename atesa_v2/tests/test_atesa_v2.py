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
from atesa_v2.configure import configure

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

    def test_init_threads_restart(self):
        """Tests successful unpickle and return of a restart.pkl file"""
        pass # todo: implement

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
        allthreads = atesa_v2.init_threads(settings)
        assert len(allthreads) == 1
        assert allthreads[0].history.init_inpcrd == ['init.rst7']
        assert allthreads[0].topology == 'topology.prmtop'

    def test_thread_get_frame_amber(self):
        """Tests thread.get_frame method with md_engine = 'amber' and frame = -1"""
        settings = configure('../../data/atesa.config')
        settings.topology = '../test_data/test.prmtop'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].jobids.append('01234')
        allthreads[0].history.prod_trajs.append(['../test_data/test.nc'])
        compare_traj = pytraj.iterload('../test_data/test.rst7', allthreads[0].topology)
        test_frame = allthreads[0].get_frame(allthreads[0].history.prod_trajs[0][0], -1, settings)
        query_traj = pytraj.iterload(test_frame, allthreads[0].topology)
        assert query_traj.n_frames == compare_traj.n_frames
        assert pytraj.center_of_mass(query_traj) == pytest.approx(pytraj.center_of_mass(compare_traj), 1e-3)

    def test_gatekeep_aimless_shooting(self):
        """Tests thread.gatekeep method with job_type = 'aimless_shooting'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.DEBUG = True
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].gatekeeper(settings) == True

    # def test_gatekeep_committor_analysis(self):
    #     """Tests thread.gatekeep method with job_type = 'committor_analysis'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'committor_analysis'
    #     settings.DEBUG = True
    #     allthreads = atesa_v2.init_threads(settings)
    #     assert allthreads[0].gatekeeper(settings) == True

    # def test_gatekeep_equilibrium_path_sampling(self):
    #     """Tests thread.gatekeep method with job_type = 'equilibrium_path_sampling'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'equilibrium_path_sampling'
    #     settings.DEBUG = True
    #     allthreads = atesa_v2.init_threads(settings)
    #     assert allthreads[0].gatekeeper(settings) == True

    def test_get_batch_file_aimless_shooting_amber(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and md_engine = 'amber'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_aimless_shooting_broken(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and invalid type = 'fwdd'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('fwdd', settings)

    # def test_get_batch_file_committor_analysis_amber(self):
    #     """Tests thread.get_batch_template with job_type = 'committor_analysis' and md_engine = 'amber'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'committor_analysis'
    #     settings.md_engine = 'amber'
    #     allthreads = atesa_v2.init_threads(settings)
    #     assert allthreads[0].get_batch_template('prod', settings) == 'amber_' + settings.batch_system + '.tpl'

    # def test_get_batch_file_committor_analysis_broken(self):
    #     """Tests thread.get_batch_template with job_type = 'committor_analysis' and invalid type = 'init'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'committor_analysis'
    #     allthreads = atesa_v2.init_threads(settings)
    #     with pytest.raises(ValueError):
    #         allthreads[0].get_batch_template('init', settings)

    # def test_get_batch_file_equilibrium_path_sampling_amber(self):
    #     """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and md_engine = 'amber'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'equilibrium_path_sampling'
    #     settings.md_engine = 'amber'
    #     allthreads = atesa_v2.init_threads(settings)
    #     assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    # def test_get_batch_file_equilibrium_path_sampling_broken(self):
    #     """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and invalid type = 'fwdd'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'equilibrium_path_sampling'
    #     allthreads = atesa_v2.init_threads(settings)
    #     with pytest.raises(ValueError):
    #         allthreads[0].get_batch_template('fwdd', settings)

    def test_get_next_step_aimless_shooting_init(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['init']
        assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    def test_get_next_step_aimless_shooting_prod(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['fwd','bwd']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['prod','prod']
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    # def test_get_next_step_committor_analysis(self):
    #     """Tests thread.get_next_step with job_type = 'committor_analysis'"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'committor_analysis'
    #     settings.committor_analysis_n = 10
    #     allthreads = atesa_v2.init_threads(settings)
    #     allthreads[0].current_type = []
    #     assert allthreads[0].get_next_step(settings) == (['prod' for null in range(settings.committor_analysis_n)], [str(int(i)) for i in range(settings.committor_analysis_n)])

    # def test_get_next_step_equilibrium_path_sampling_init(self):
    #     """Tests thread.get_next_step with job_type = 'equilibrium_path_sampling' and type = ['init']"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'equilibrium_path_sampling'
    #     allthreads = atesa_v2.init_threads(settings)
    #     allthreads[0].current_type = ['init']
    #     assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    # def test_get_next_step_equilibrium_path_sampling_prod(self):
    #     """Tests thread.get_next_step with job_type = 'equilibrium_path_sampling' and type = ['fwd','bwd']"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'equilibrium_path_sampling'
    #     allthreads = atesa_v2.init_threads(settings)
    #     allthreads[0].current_type = ['prod','prod']
    #     assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
