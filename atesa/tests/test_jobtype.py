"""
Unit and regression test for jobtype.py.
"""

# Import package, test suite, and other packages as needed
# import atesa
import pytest
import sys
import glob
import os
import filecmp
import shutil
from atesa.configure import configure
from atesa import factory
from atesa import main

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa/tests/test_temp'):
                os.mkdir('atesa/tests/test_temp')
            os.chdir('atesa/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_check_termination_aimless_shooting_init(self):
        """Tests check_termination with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False  # never terminate after an 'init' step

    def test_check_termination_aimless_shooting_prod(self):
        """Tests check_termination with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False  # todo: update after implementation
        assert allthreads[0].terminated == False              # todo: update after implementation


    def test_update_results_aimless_shooting_init(self):
        """Tests update_results with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_results(allthreads[0], allthreads, settings)
        assert os.path.exists('restart.pkl')

    def test_update_results_aimless_shooting_prod(self):
        """Tests update_results with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.5, 2.0], ['lt', 'gt']]
        settings.commit_bwd = [[1, 2], [3, 4], [2.0, 1.5], ['gt', 'lt']]
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[0].history.init_coords = [['test_velocities.rst7_0_init.rst7']]
        shutil.copy('../test_data/test_velocities.rst7', 'test_velocities.rst7_0_init.rst7')  # create the needed file
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_results(allthreads[0], allthreads, settings)
        assert os.path.exists('restart.pkl')
        assert os.path.exists('status.txt')
        assert allthreads[0].history.prod_results == [['', '']]

    def test_algorithm_aimless_shooting_init(self):
        """Tests algorithm with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        allthreads[0].history.init_coords = [['test_velocities.rst7_0_init.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert allthreads[0].current_type == []     # result for missing .rst7 file (haven't copied it yet)
        allthreads[0].current_type = ['init']       # reset last result
        shutil.copy('../test_data/test_velocities.rst7', 'test_velocities.rst7_0_init.rst7')    # create the needed file
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert allthreads[0].current_type == ['init']   # results for .rst7 was found
        assert allthreads[0].history.init_coords == [['test_velocities.rst7_0_init.rst7', 'test_velocities.rst7_0_init_bwd.rst7']]

    def test_algorithm_aimless_shooting_prod_not_always_new_not_accepted(self):
        """Tests algorithm with job_type = 'aimless_shooting', always_new = False and thread.current_type =
        ['prod', 'prod'] for a move that isn't accepted"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.always_new = False
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results.append(['fwd', 'fwd'])  # not an accepted move
        allthreads[0].history.prod_trajs.append(['../test_data/test.nc', '../test_data/test.nc'])
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert allthreads[0].history.init_inpcrd[-1] == allthreads[0].history.init_inpcrd[-2]

    def test_algorithm_aimless_shooting_prod_always_new_not_accepted(self):
        """Tests algorithm with job_type = 'aimless_shooting', always_new = True and thread.current_type =
        ['prod', 'prod'] for a move that isn't accepted"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.topology = '../test_data/test.prmtop'
        settings.always_new = True
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        settings.min_dt = -1
        settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results = [['bwd', 'fwd'], ['fwd', 'fwd']]      # accepted then not accepted
        allthreads[0].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc'], ['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[0].suffix = 1
        allthreads[0].history.last_accepted = 0
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert filecmp.cmp(allthreads[0].history.init_inpcrd[1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc
        os.remove('../test_data/test.nc_frame_-1.rst7') # have to do this manually because aimless shooting's mdengine getframe method keeps the '../test_data/' in front of the file name

    def test_algorithm_aimless_shooting_prod_accepted(self):
        """Tests algorithm with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod'] for an accepted move"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.topology = '../test_data/test.prmtop'
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        settings.degeneracy = 1
        settings.min_dt = -1
        settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results = [['bwd', 'bwd'], ['fwd', 'bwd']]  # not accepted then accepted
        allthreads[0].history.prod_trajs = [['not_a_real_file.nc', 'not_a_real_file.nc'], ['../test_data/test.nc', '../test_data/test.nc']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert filecmp.cmp(allthreads[0].history.init_inpcrd[1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc
        os.remove('../test_data/test.nc_frame_-1.rst7') # have to do this manually because aimless shooting's mdengine getframe method keeps the '../test_data/' in front of the file name

    def test_update_history_aimless_shooting_init(self):
        """Tests update_history with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        these_kwargs = {'rst': 'fakey_mcfakename.rst'}
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.init_coords[-1] == ['fakey_mcfakename.rst']

    def test_update_history_aimless_shooting_prod(self):
        """Tests update_history with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        these_kwargs = {'nc': 'fakey_mcfakename.nc'}
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs[-1] == ['fakey_mcfakename.nc']
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs[-1] == ['fakey_mcfakename.nc', 'fakey_mcfakename.nc']

    def test_get_inpcrd_aimless_shooting_init(self):
        """Tests get_inpcrd with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_inpcrd(allthreads[0]) == allthreads[0].history.init_inpcrd
        assert allthreads[0].history.init_inpcrd == ['test_velocities.rst7']

    def test_get_inpcrd_aimless_shooting_prod(self):
        """Tests get_inpcrd with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.init_coords = [['not_a_real_file_at_all_init.rst7', 'not_a_real_file_at_all_init_bwd.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_inpcrd(allthreads[0]) == ['not_a_real_file_at_all_init.rst7', 'not_a_real_file_at_all_init_bwd.rst7']

    def test_get_initial_coordinates_aimless_shooting(self):
        """Tests get_initial_coordinates with job_type = 'aimless_shooting'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.initial_coordinates = ['../test_data/test.rst7', '../test_data/test_two_init.rst7']
        settings.degeneracy = 2
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_initial_coordinates(settings) == ['test.rst7_0', 'test.rst7_1', 'test_two_init.rst7_0', 'test_two_init.rst7_1']

    def test_gatekeep_aimless_shooting(self):
        """Tests gatekeep with job_type = 'aimless_shooting'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.DEBUG = True
        allthreads = main.init_threads(settings)
        allthreads[0].jobids = ['123456']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.gatekeeper(allthreads[0], settings) == True

    def test_check_for_successful_step_aimless_shooting_init(self):
        """Tests check_for_successful_step with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        allthreads[0].history.init_coords = [['some_init_coords.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_for_successful_step(allthreads[0], settings) == False    # necessary file does not yet exist
        shutil.copy('../test_data/test.rst7', 'some_init_coords.rst7')      # make the necessary file
        assert jobtype.check_for_successful_step(allthreads[0], settings) == True     # necessary file exists

    def test_get_next_step_aimless_shooting_init(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    def test_get_next_step_aimless_shooting_prod(self):
        """Tests thread.get_next_step with job_type = 'aimless_shooting' and type = ['fwd','bwd']"""
        settings = configure('../../data/atesa.config')
        settings.degeneracy = 1
        settings.job_type = 'aimless_shooting'
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod','prod']
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    def test_check_for_successful_step_aimless_shooting_prod(self):
        """Tests check_for_successful_step with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_trajs = [['some_prod_traj_1.nc', 'some_prod_traj_2.nc']]
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_for_successful_step(allthreads[0], settings) == False    # necessary files do not yet exist
        shutil.copy('../test_data/test.nc', 'some_prod_traj_1.nc')          # make one necessary file
        assert jobtype.check_for_successful_step(allthreads[0], settings) == False    # still missing one
        shutil.copy('../test_data/test.nc', 'some_prod_traj_2.nc')          # make other necessary file
        assert jobtype.check_for_successful_step(allthreads[0], settings) == True     # both files exist

    def test_get_batch_file_aimless_shooting_amber(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and md_engine = 'amber'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.md_engine = 'amber'
        allthreads = main.init_threads(settings)
        assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_aimless_shooting_broken(self):
        """Tests thread.get_batch_template with job_type = 'aimless_shooting' and invalid type = 'fwdd'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        allthreads = main.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('fwdd', settings)

    def test_check_termination_committor_analysis(self):
        """Tests check_termination with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod', 'prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False
        assert allthreads[0].terminated == True

    def test_update_results_committor_analysis(self):
        """Tests update_results with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        settings.committor_analysis_use_rc_out = False
        settings.topology = '../test_data/test.prmtop'
        settings.commit_fwd = [[1, 2], [3, 4], [1.5, 2.0], ['lt', 'gt']]
        settings.commit_bwd = [[1, 2], [3, 4], [2.0, 1.5], ['gt', 'lt']]
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod', 'prod']
        allthreads[0].history.prod_trajs = ['../test_data/test.nc', '../test_data/test.nc', '../test_data/test.nc']
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_results(allthreads[0], allthreads, settings)
        assert os.path.exists('committor_analysis.out')
        assert allthreads[0].history.prod_results == ['', '', '']

    def test_update_history_committor_analysis(self):
        """Tests update_history with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        settings.committor_analysis_use_rc_out = False
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod', 'prod']
        these_kwargs = {'nc': 'fakey_mcfakename.nc'}
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs == ['fakey_mcfakename.nc']
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs == ['fakey_mcfakename.nc', 'fakey_mcfakename.nc']

    def test_get_inpcrd_committor_analysis(self):
        """Tests get_inpcrd with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        settings.committor_analysis_use_rc_out = False
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod', 'prod']
        allthreads[0].history.prod_inpcrd = ['not_a_real_file_at_all_init.rst7']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_inpcrd(allthreads[0]) == ['not_a_real_file_at_all_init.rst7', 'not_a_real_file_at_all_init.rst7', 'not_a_real_file_at_all_init.rst7']

    def test_get_initial_coordinates_committor_analysis_rc_out(self):
        """Tests get_initial_coordinates with job_type = 'committor_analysis' using an RC out file"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = True
        settings.path_to_rc_out = '../test_data/rc_eval.out'
        settings.rc_threshold = 0.002
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_initial_coordinates(settings) == ['1.1_1.2_2.2_2.6.rst7_2_10_init_fwd.rst_2_1_2.rst7_5_init_fwd.rst', '1.1_1.2_2.2_2.6.rst7_2_10_init_fwd.rst_4_58_3.rst7_4_init_fwd.rst']

    def test_get_initial_coordinates_committor_analysis_rc_out_does_not_exist(self):
        """Tests get_initial_coordinates with job_type = 'committor_analysis' using an RC out file that does not exist"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = True
        settings.path_to_rc_out = '../test_data/rc_eval_FAKE.out'
        settings.rc_threshold = 0.002
        jobtype = factory.jobtype_factory(settings.job_type)
        with pytest.raises(FileNotFoundError):
            allthreads = main.init_threads(settings)

    def test_get_initial_coordinates_committor_analysis_rc_out_no_shooting_points(self):
        """Tests get_initial_coordinates with job_type = 'committor_analysis' using an RC out file containing no
        shooting points within the threshold"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = True
        settings.path_to_rc_out = '../test_data/rc_eval.out'
        settings.rc_threshold = 0.0002
        jobtype = factory.jobtype_factory(settings.job_type)
        with pytest.raises(RuntimeError):
            allthreads = main.init_threads(settings)

    def test_gatekeep_committor_analysis(self):
        """Tests gatekeep with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.DEBUG = True
        settings.committor_analysis_use_rc_out = False
        allthreads = main.init_threads(settings)
        allthreads[0].jobids = ['123456']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.gatekeeper(allthreads[0], settings) == True

    def test_check_for_successful_step_committor_analysis(self):
        """Tests check_for_successful_step with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        allthreads = main.init_threads(settings)
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_for_successful_step(allthreads[0], settings) == True

    def test_get_next_step_committor_analysis(self):
        """Tests get_next_step with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        settings.committor_analysis_n = 3
        allthreads = main.init_threads(settings)
        jobtype = factory.jobtype_factory(settings.job_type)
        this_types, this_names = jobtype.get_next_step(allthreads[0], settings)
        assert this_types == ['prod', 'prod', 'prod']
        assert this_names == ['0', '1', '2']

    def test_get_batch_template_committor_analysis(self):
        """Tests get_batch_template with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.md_engine = 'amber'
        settings.batch_system = 'slurm'
        settings.committor_analysis_use_rc_out = False
        settings.path_to_templates = sys.path[0] + '/atesa/data/templates'
        allthreads = main.init_threads(settings)
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_batch_template(allthreads[0], 'prod', settings) == 'amber_slurm.tpl'

    def test_check_termination_equilibrium_path_sampling_global(self):
        """Tests check_termination global criterion with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        settings.samples_per_window = -1
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False
        assert allthreads[0].terminated == False

    def test_check_termination_equilibrium_path_sampling_thread(self):
        """Tests check_termination thread criterion with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[0].history.prod_lens = [[2, 2]]
        allthreads[0].history.init_coords = [['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_results(allthreads[0], allthreads, settings)
        settings.samples_per_window = 1
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False  # global termination False
        assert allthreads[0].terminated == True     # thread termination True

    def test_update_results_equilibrium_path_sampling(self):
        """Tests update_results with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        settings.initial_coordinates = ['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_trajs = [['../test_data/test.nc', '../test_data/test.nc']]
        allthreads[0].history.prod_lens = [[2, 2]]
        allthreads[0].history.init_coords = [['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_results(allthreads[0], allthreads, settings)
        assert os.path.exists('eps.out')
        assert allthreads[0].history.prod_results[-1] == pytest.approx([-34.29, -35.89, -35.17, -35.89, -35.17], 1E-3)

    def test_update_history_equilibrium_path_sampling_prod(self):
        """Tests update_history with job_type = 'equilibrium_path_sampling' with current_type = ['prod', 'prod']"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        these_kwargs = {'nc': 'fakey_mcfakename.nc'}
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs == [['fakey_mcfakename.nc']]
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.prod_trajs == [['fakey_mcfakename.nc', 'fakey_mcfakename.nc']]

    def test_update_history_equilibrium_path_sampling_init(self):
        """Tests update_history with job_type = 'equilibrium_path_sampling' with current_type = ['init']"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        these_kwargs = {'rst': 'fakey_mcfakename.rst7'}
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.update_history(allthreads[0], settings, **these_kwargs)
        assert allthreads[0].history.init_coords == [['fakey_mcfakename.rst7']]

    def test_update_history_equilibrium_path_sampling_out_of_bounds(self):
        """Tests update_history with job_type = 'equilibrium_path_sampling' with an initially out-of-bounds RC value"""
        shutil.copy('../../data/atesa.config', 'eps.config')
        with open('eps.config', 'a') as f:  # need to set these things before calling configure()
            f.write('\njob_type = \'equilibrium_path_sampling\'')
            f.write('\neps_rc_min = -0.1')  # crazy wide range so everything gets included
            f.write('\neps_rc_max = 0.1')
            f.write('\neps_rc_step = 0.1')
            f.write('\neps_overlap = 0.01')
            f.write('\neps_dynamic_seed = 3')
            f.close()
        settings = configure('eps.config')
        settings.DEBUG = True
        settings.job_type = 'equilibrium_path_sampling'
        settings.topology = '../test_data/test.prmtop'
        settings.rc_reduced_cvs = False
        settings.include_qdot = False
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
        settings.rc_definition = '1.00 + 2.34*CV1 - 0.67*CV2'
        settings.initial_coordinates = ['../test_data/test.rst7']
        with pytest.raises(RuntimeError):
            allthreads = main.init_threads(settings)

    def test_get_inpcrd_equilibrium_path_sampling(self):
        """Tests get_inpcrd with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.init_coords = [['not_a_real_file_at_all_init.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_inpcrd(allthreads[0]) == ['not_a_real_file_at_all_init.rst7']

    def test_get_initial_coordinates_equilibrium_path_sampling(self):
        """Tests get_initial_coordinates with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_initial_coordinates(settings) == ['test.rst7']

    def test_gatekeep_equilibrium_path_sampling(self):
        """Tests gatekeep with job_type = 'committor_analysis'"""
        settings = config_equilibrium_path_sampling()
        settings.DEBUG = True
        allthreads = main.init_threads(settings)
        allthreads[0].jobids = ['123456']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.gatekeeper(allthreads[0], settings) == True

    def test_check_for_successful_step_equilibrium_path_sampling(self):
        """Tests check_for_successful_step with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        allthreads[0].history.init_coords = [['test_velocities.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_for_successful_step(allthreads[0], settings) == False
        shutil.copy('../test_data/test_velocities.rst7', 'test_velocities.rst7')
        assert jobtype.check_for_successful_step(allthreads[0], settings) == True
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_trajs = [['test.nc', 'test.nc']]
        assert jobtype.check_for_successful_step(allthreads[0], settings) == False
        shutil.copy('../test_data/test.nc', 'test.nc')
        assert jobtype.check_for_successful_step(allthreads[0], settings) == True

    def test_get_next_step_equilibrium_path_sampling_first(self):
        """Tests get_next_step with job_type = 'equilibrium_path_sampling' and an empty current_type"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = []
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    def test_get_next_step_equilibrium_path_sampling_init(self):
        """Tests get_next_step with job_type = 'equilibrium_path_sampling' and current_type = ['init']"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        assert allthreads[0].get_next_step(settings) == (['prod', 'prod'], ['fwd', 'bwd'])

    def test_get_next_step_equilibrium_path_sampling_prod(self):
        """Tests get_next_step with job_type = 'equilibrium_path_sampling' and current_type = ['prod', 'prod']"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        assert allthreads[0].get_next_step(settings) == (['init'], ['init'])

    def test_get_batch_file_equilibrium_path_sampling_amber(self):
        """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and md_engine = 'amber'"""
        settings = config_equilibrium_path_sampling()
        settings.job_type = 'equilibrium_path_sampling'
        settings.md_engine = 'amber'
        allthreads = main.init_threads(settings)
        assert allthreads[0].get_batch_template('init', settings) == 'amber_' + settings.batch_system + '.tpl'

    def test_get_batch_file_equilibrium_path_sampling_broken(self):
        """Tests thread.get_batch_template with job_type = 'equilibrium_path_sampling' and invalid type = 'fwdd'"""
        settings = config_equilibrium_path_sampling()
        settings.job_type = 'equilibrium_path_sampling'
        allthreads = main.init_threads(settings)
        with pytest.raises(ValueError):
            allthreads[0].get_batch_template('fwdd', settings)

    def test_algorithm_equilibrium_path_sampling_init(self):
        """Tests algorithm with job_type = 'equilibrium_path_sampling' and thread.current_type = ['init']"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        allthreads[0].history.init_coords = [['test_velocities.rst7_0_init.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert allthreads[0].current_type == []     # result for missing .rst7 file (haven't copied it yet)
        allthreads[0].current_type = ['init']       # reset last result
        shutil.copy('../test_data/test_velocities.rst7', 'test_velocities.rst7_0_init.rst7')    # create the needed file
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert allthreads[0].current_type == ['init']   # results for .rst7 was found
        assert allthreads[0].history.init_coords == [['test_velocities.rst7_0_init.rst7', 'test_velocities.rst7_0_init_bwd.rst7']]

    def test_algorithm_equilibrium_path_sampling_prod_not_accepted(self):
        """Tests algorithm with job_type = 'equilibrium_path_sampling' and thread.current_type = ['prod', 'prod'] for a
        move that isn't accepted"""
        settings = config_equilibrium_path_sampling()
        settings.min_dt = -1
        settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
        settings.eps_dynamic_seed = 0
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results = [[-34.12, -35.1, -36], [12, 13, 11]]      # accepted then not accepted
        allthreads[0].history.prod_trajs = [['test.nc', 'test.nc'], ['not_a_real_file.nc', 'not_a_real_file.nc']]
        allthreads[0].history.init_coords = [['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']]
        allthreads[0].suffix = 1
        allthreads[0].history.last_accepted = 0
        shutil.copy('../test_data/test.nc', 'test.nc')
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert filecmp.cmp(allthreads[0].history.init_inpcrd[1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc

    def test_algorithm_equilibrium_path_sampling_prod_not_accepted_no_accepted_yet(self):
        """Tests algorithm with job_type = 'equilibrium_path_sampling' and thread.current_type = ['prod', 'prod'] for a
        move that isn't accepted and with no accepted moves in the thread's history"""
        settings = config_equilibrium_path_sampling()
        settings.min_dt = -1
        settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results = [[16, 13, 15], [12, 13, 11]]      # accepted then not accepted
        allthreads[0].history.prod_trajs = [['test.nc', 'test.nc'], ['test.nc', 'test.nc']]
        allthreads[0].history.init_coords = [['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']]
        allthreads[0].suffix = 1
        shutil.copy('../test_data/test.nc', 'test.nc')
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert filecmp.cmp(allthreads[0].history.init_inpcrd[1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc

    def test_algorithm_equilibrium_path_sampling_prod_accepted(self):
        """Tests algorithm with job_type = 'equilibrium_path_sampling' and thread.current_type = ['prod', 'prod'] for
        an accepted move"""
        settings = config_equilibrium_path_sampling()
        settings.min_dt = -1
        settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results = [[12, 13, 11], [-34.12, -35.1, -36]]  # not accepted then accepted
        allthreads[0].history.prod_trajs = [['not_a_real_file.nc', 'not_a_real_file.nc'], ['test.nc', 'test.nc']]
        allthreads[0].history.init_coords = [['../test_data/test_velocities.rst7', '../test_data/test_velocities.rst7']]
        shutil.copy('../test_data/test.nc', 'test.nc')
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, allthreads, settings)
        assert filecmp.cmp(allthreads[0].history.init_inpcrd[1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc

    def test_get_input_file_aimless_shooting(self):
        """Tests get_input_file with job_type = 'aimless_shooting'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == settings.path_to_input_files + '/' + settings.job_type + '_' + allthreads[0].current_type[0] + '_' + settings.md_engine + '.in'

    def test_get_input_file_committor_analysis(self):
        """Tests get_input_file with job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == settings.path_to_input_files + '/' + settings.job_type + '_' + allthreads[0].current_type[0] + '_' + settings.md_engine + '.in'

    def test_get_input_file_equilibrium_path_sampling(self):
        """Tests get_input_file with job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == settings.path_to_input_files + '/' + settings.job_type + '_' + allthreads[0].current_type[0] + '_' + settings.md_engine + '.in'
        allthreads[0].current_type = ['prod', 'prod']
        new_file = jobtype.get_input_file(allthreads[0], 0, settings)
        assert os.path.exists(new_file)
        assert allthreads[0].history.prod_lens[-1][0] + allthreads[0].history.prod_lens[-1][1] == settings.eps_n_steps
        new_file = jobtype.get_input_file(allthreads[0], 1, settings)
        assert os.path.exists(new_file)

    def test_get_input_file_rxncor_umbrella_sampling(self):
        """Tests get_input_file with job_type = 'umbrella_sampling'"""
        settings = configure('../../data/atesa.config')
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', '(mdtraj.compute_distances(mtraj, numpy.array([[7177, 7178]]))[0][0] * 10) - (mdtraj.compute_distances(mtraj, numpy.array([[4272, 7178]]))[0][0] * 10)']
        settings.us_rc_min = -2
        settings.us_rc_max = 2
        settings.us_rc_step = 0.5
        settings.us_degeneracy = 2
        settings.us_implementation = 'amber_rxncor'
        settings.job_type = 'umbrella_sampling'
        settings.as_out_file = '../test_data/as.out'
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['prod']
        allthreads[0].history.window = 2.5
        allthreads[0].history.index = 1
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == 'umbrella_sampling_2.5_1.in'
        assert filecmp.cmp('umbrella_sampling_2.5_1.in', '../test_data/umbrella_sampling_2.5_1.in')

    @classmethod
    def teardown_method(self, method):
        "Runs at end of each method"
        for filename in glob.glob(sys.path[0] + '/atesa/tests/test_temp/*'):
            os.remove(filename)

def config_equilibrium_path_sampling():
    """Sets up configuration settings for equilibrium path sampling tests, to be overwritten as needed"""
    shutil.copy('../../data/atesa.config', 'eps.config')
    with open('eps.config', 'a') as f:  # need to set these things before calling configure()
        f.write('\njob_type = \'equilibrium_path_sampling\'')
        f.write('\neps_rc_min = -50')   # crazy wide range so everything gets included
        f.write('\neps_rc_max = 50')
        f.write('\neps_rc_step = 1')
        f.write('\neps_rc_overlap = 0.1')
        f.write('\neps_n_steps = 6')
        f.write('\neps_out_freq = 1')
        f.write('\neps_dynamic_seed = 3')
        f.write('\nsamples_per_window = -1')
        f.close()

    settings = configure('eps.config')
    settings.DEBUG = True
    settings.job_type = 'equilibrium_path_sampling'
    settings.topology = '../test_data/test.prmtop'
    settings.rc_reduced_cvs = False
    settings.include_qdot = False
    settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', 'pytraj.angle(traj, \'@2 @3 @4\')[0]']
    settings.rc_definition = '1.00 + 2.34*CV1 - 0.67*CV2'
    settings.initial_coordinates = ['../test_data/test.rst7']

    return settings
