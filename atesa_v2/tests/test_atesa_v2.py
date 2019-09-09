"""
Unit and regression test for the atesa_v2 package.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import pytest
import sys
from ..configure import configure
import pytraj
import os
import glob

class Tests(object):
    def test_atesa_v2_imported(self):
        """Sample test, will always pass so long as import statement worked"""
        assert "atesa_v2" in sys.modules

    def test_init_threads_restart(self):
        """Tests successful unpickle and return of a restart.pkl file"""
        pass # todo: implement

    def test_configure_broken(self):
        """Tests configure.py with a non-existent file"""
        with pytest.raises(FileNotFoundError):
            settings = configure('atesa_v2/data/ates.config')

    def test_init_threads_new(self):
        """Tests successful initialization of new threads"""
        settings = configure('atesa_v2/data/atesa.config')
        allthreads = atesa_v2.init_threads(settings)
        assert len(allthreads) == 1
        assert allthreads[0].coordinates == 'init.rst7'
        assert allthreads[0].topology == 'topology.prmtop'

    def test_thread_get_last_frame_amber(self):
        """Tests thread.get_last_frame method with md_engine = 'amber'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.topology = 'atesa_v2/tests/test_data/test.prmtop'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].jobids.append('01234')
        allthreads[0].traj_files.append('atesa_v2/tests/test_data/test.nc')
        compare_traj = pytraj.iterload('atesa_v2/tests/test_data/test.rst7', allthreads[0].topology)
        query_traj = pytraj.iterload(allthreads[0].get_last_frame(0, settings), allthreads[0].topology)
        assert query_traj.n_frames == compare_traj.n_frames
        assert pytraj.center_of_mass(query_traj) == pytest.approx(pytraj.center_of_mass(compare_traj), 1e-3)
        os.remove('atesa_v2/tests/test_data/test.nc_last_frame.rst7')

    def test_check_commit(self):
        """Tests check_commit using a dummy coordinate file"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.topology = 'atesa_v2/tests/test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == 'fwd'
        settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.7], ['gt', 'lt']]
        settings.commit_fwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == 'bwd'
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [0.5, 2.0], ['lt', 'gt']]
        assert atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings) == ''
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
        with pytest.raises(ValueError):
            atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings)
        settings.commit_fwd = [[1, 2], [3, 4], [1.0, 1.5], ['gt', 'lt']]
        settings.commit_bwd = [[1, 2], [3, 4], [1.0, 1.5], ['t', 'lt']]
        with pytest.raises(ValueError):
            atesa_v2.check_commit('atesa_v2/tests/test_data/test.rst7', settings)

    def test_gatekeep_aimless_shooting(self):
        """Tests thread.gatekeep method with job_type = 'aimless_shooting'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.DEBUG = True
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].gatekeeper(settings) == True

    def test_gatekeep_committor_analysis(self):
        """Tests thread.gatekeep method with job_type = 'committor_analysis'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.DEBUG = True
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].gatekeeper(settings) == True

    def test_gatekeep_equilibrium_path_sampling(self):
        """Tests thread.gatekeep method with job_type = 'equilibrium_path_sampling'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'equilibrium_path_sampling'
        settings.DEBUG = True
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].gatekeeper(settings) == True

    def test_get_batch_file_aimless_shooting_amber(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and md_engine = 'amber'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_aimless_shooting_broken(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and invalid type = 'fwdd'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('fwdd', settings)

    def test_get_batch_file_committor_analysis_amber(self):
        """Tests thread.get_batch_template with job_type = 'committor_analysis' and md_engine = 'amber'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].get_batch_template('prod', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_committor_analysis_broken(self):
        """Tests thread.get_batch_template with job_type = 'committor_analysis' and invalid type = 'init'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'committor_analysis'
        allthreads = atesa_v2.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('init', settings)

    def test_get_batch_file_equilibrium_path_sampling_amber(self):
        """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and md_engine = 'amber'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'equilibrium_path_sampling'
        settings.md_engine = 'amber'
        allthreads = atesa_v2.init_threads(settings)
        assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_equilibrium_path_sampling_broken(self):
        """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and invalid type = 'fwdd'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'equilibrium_path_sampling'
        allthreads = atesa_v2.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('fwdd', settings)

    def test_get_next_step_aimless_shooting_init(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['init']"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['init']
        assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    def test_get_next_step_aimless_shooting_prod(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['fwd','bwd']"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['prod','prod']
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    def test_get_next_step_committor_analysis(self):
        """Tests thread.get_next_step with job_type = 'committor_analysis'"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_n = 10
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = []
        assert allthreads[0].get_next_step(settings) == (['prod' for null in range(settings.committor_analysis_n)], [str(int(i)) for i in range(settings.committor_analysis_n)])

    def test_get_next_step_equilibrium_path_sampling_init(self):
        """Tests thread.get_next_step with job_type = 'equilibrium_path_sampling' and type = ['init']"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'equilibrium_path_sampling'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['init']
        assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    def test_get_next_step_equilibrium_path_sampling_prod(self):
        """Tests thread.get_next_step with job_type = 'equilibrium_path_sampling' and type = ['fwd','bwd']"""
        settings = configure('atesa_v2/data/atesa.config')
        settings.job_type = 'equilibrium_path_sampling'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['prod','prod']
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    @classmethod
    def teardown_class(cls):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
