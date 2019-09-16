"""
Unit and regression test for jobtype.py.
"""

# Import package, test suite, and other packages as needed
import atesa_v2
import pytest
from atesa_v2.configure import configure
try:
    import factory
except ModuleNotFoundError:
    import atesa_v2.factory as factory
import sys
import glob
import os
import filecmp
import shutil

class Tests(object):
    def setup_method(self, test_method):
        try:
            if not os.path.exists('atesa_v2/tests/test_temp'):
                os.mkdir('atesa_v2/tests/test_temp')
            os.chdir('atesa_v2/tests/test_temp')
        except FileNotFoundError:
            pass

    def test_check_termination_aimless_shooting_init(self):
        """Tests check_termination with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False  # never terminate after an 'init' step

    def test_check_termination_aimless_shooting_prod(self):
        """Tests check_termination with job_type = 'aimless_shooting' and thread.current_type = ['prod', 'prod']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.check_termination(allthreads[0], allthreads, settings) == False  # todo: update after implementation
        assert allthreads[0].terminated == False                                        # todo: update after implementation

    def test_update_results_aimless_shooting_init(self):
        """Tests update_results with job_type = 'aimless_shooting' and thread.current_type = ['init']"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        allthreads = atesa_v2.init_threads(settings)
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
        allthreads = atesa_v2.init_threads(settings)
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
        settings.initial_coordinates = ['test_velocities.rst7']
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['init']
        allthreads[0].history.init_coords = [['test_velocities.rst7_0_init.rst7']]
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, settings)
        assert allthreads[0].current_type == []     # result for missing .rst7 file (haven't copied it yet)
        allthreads[0].current_type = ['init']       # reset last result
        shutil.copy('../test_data/test_velocities.rst7', 'test_velocities.rst7_0_init.rst7')    # create the needed file
        jobtype.algorithm(allthreads[0], allthreads, settings)
        assert allthreads[0].current_type == ['init']   # results for .rst7 was found
        assert allthreads[0].history.init_coords == [['test_velocities.rst7_0_init.rst7', 'test_velocities.rst7_0_init_bwd.rst7']]

    def test_algorithm_aimless_shooting_prod_not_always_new_not_accepted(self):
        """Tests algorithm with job_type = 'aimless_shooting', always_new = False and thread.current_type =
        ['prod', 'prod'] for a move that isn't accepted"""
        settings = configure('../../data/atesa.config')
        settings.job_type = 'aimless_shooting'
        settings.always_new = False
        settings.initial_coordinates = ['test_velocities.rst7']
        allthreads = atesa_v2.init_threads(settings)
        allthreads[0].current_type = ['prod', 'prod']
        allthreads[0].history.prod_results.append(['fwd', 'fwd'])      # not an accepted move
        jobtype = factory.jobtype_factory(settings.job_type)
        jobtype.algorithm(allthreads[0], allthreads, settings)
        assert allthreads[0].history.init_inpcrd[-1] == allthreads[0].history.init_inpcrd[-2]

    # def test_algorithm_aimless_shooting_prod_always_new_not_accepted(self):
    #     """Tests algorithm with job_type = 'aimless_shooting', always_new = True and thread.current_type =
    #     ['prod', 'prod'] for a move that isn't accepted"""
    #     settings = configure('../../data/atesa.config')
    #     settings.job_type = 'aimless_shooting'
    #     settings.always_new = False
    #     settings.initial_coordinates = ['test_velocities.rst7']
    #     settings.min_dt = -1
    #     settings.max_dt = -1     # set these to the same value to guarantee which frame is chosen
    #     allthreads = atesa_v2.init_threads(settings)
    #     allthreads[0].current_type = ['prod', 'prod']
    #     allthreads[0].history.prod_results = [['bwd', 'fwd'], ['fwd', 'fwd']]      # accepted then not accepted
    #     allthreads[0].history.prod_trajs = [['test.nc', 'test.nc'], ['not_a_real_file.nc', 'not_a_real_file.nc']]
    #     jobtype = factory.jobtype_factory(settings.job_type)
    #     jobtype.algorithm(allthreads[0], allthreads, settings)
    #     assert filecmp.cmp(allthreads[0].history.init_inpcrd[-1], '../test_data/test.rst7') # test.rst7 is last frame of test.nc

    @classmethod
    def teardown_method(self, method):
        "Runs at end of class"
        for filename in glob.glob(sys.path[0] + '/atesa_v2/tests/test_temp/*'):
            os.remove(filename)
