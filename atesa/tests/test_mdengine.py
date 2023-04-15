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
import numpy
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

    def test_get_input_file_amber_aimless_shooting(self):
        """Tests get_input_file with md_engine = 'amber' and job_type = 'aimless_shooting'"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'amber'
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == settings.path_to_input_files + '/' + settings.job_type + '_' + allthreads[0].current_type[0] + '_' + settings.md_engine + '.in'

    def test_get_input_file_amber_committor_analysis(self):
        """Tests get_input_file with md_engine = 'amber' and job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'amber'
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        allthreads[0].current_type = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        assert jobtype.get_input_file(allthreads[0], 0, settings) == settings.path_to_input_files + '/' + settings.job_type + '_' + allthreads[0].current_type[0] + '_' + settings.md_engine + '.in'

    def test_get_input_file_amber_equilibrium_path_sampling(self):
        """Tests get_input_file with md_engine = 'amber' and job_type = 'equilibrium_path_sampling'"""
        settings = config_equilibrium_path_sampling()
        settings.md_engine = 'amber'
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

    def test_get_input_file_amber_rxncor_umbrella_sampling(self):
        """Tests get_input_file with md_engine = 'amber' and job_type = 'umbrella_sampling'"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'amber'
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

    def test_get_input_file_cp2k_aimless_shooting_init(self):
        """Tests get_input_file with md_engine = 'cp2k' and job_type = 'aimless_shooting' for an init step"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'cp2k'
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        thread = allthreads[0]
        thread.current_type = ['init']
        thread.current_name = ['init']
        jobtype = factory.jobtype_factory(settings.job_type)
        job_index = 0
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]
        this_inpcrd = ['../test_data/test_velocities.rst7']
        inp_kwargs = {'name': thread.name + '_' + name,
                         'nodes': eval('settings.' + type + '_nodes'),
                         'taskspernode': eval('settings.' + type + '_ppn'),
                         'walltime': eval('settings.' + type + '_walltime'),
                         'mem': eval('settings.' + type + '_mem'),
                         'solver': eval('settings.' + type + '_solver'),
                         'out': thread.name + '_' + name + '.out',
                         'prmtop': thread.topology,
                         'inpcrd': this_inpcrd[job_index],
                         'rst': thread.name + '_' + name + '.rst7',
                         'nc': thread.name + '_' + name + '.nc',
                         'working_directory': settings.working_directory,
                         'extra': eval('settings.' + type + '_extra')}
        numpy.random.seed(0)
        inp = jobtype.get_input_file(thread, job_index, settings, **inp_kwargs)
        assert filecmp.cmp(inp, '../test_data/cp2k_as_test_velocities.rst7_0_init.inp')

    def test_get_input_file_cp2k_aimless_shooting_prod(self):
        """Tests get_input_file with md_engine = 'cp2k' and job_type = 'aimless_shooting' for a prod step"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'cp2k'
        settings.job_type = 'aimless_shooting'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        thread = allthreads[0]
        thread.current_type = ['prod', 'prod']
        thread.current_name = ['fwd', 'bwd']
        jobtype = factory.jobtype_factory(settings.job_type)
        job_index = 0
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]
        this_inpcrd = ['../test_data/test_velocities.rst7']
        inp_kwargs = {'name': thread.name + '_' + name,
                         'nodes': eval('settings.' + type + '_nodes'),
                         'taskspernode': eval('settings.' + type + '_ppn'),
                         'walltime': eval('settings.' + type + '_walltime'),
                         'mem': eval('settings.' + type + '_mem'),
                         'solver': eval('settings.' + type + '_solver'),
                         'out': thread.name + '_' + name + '.out',
                         'prmtop': thread.topology,
                         'inpcrd': this_inpcrd[job_index],
                         'rst': thread.name + '_' + name + '.rst7',
                         'nc': thread.name + '_' + name + '.nc',
                         'working_directory': settings.working_directory,
                         'extra': eval('settings.' + type + '_extra')}
        numpy.random.seed(0)
        inp = jobtype.get_input_file(thread, job_index, settings, **inp_kwargs)
        assert filecmp.cmp(inp, '../test_data/cp2k_as_test_velocities.rst7_0_fwd.inp')

    def test_get_input_file_cp2k_committor_analysis(self):
        """Tests get_input_file with md_engine = 'cp2k' and job_type = 'committor_analysis'"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'cp2k'
        settings.job_type = 'committor_analysis'
        settings.committor_analysis_use_rc_out = False
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        allthreads = main.init_threads(settings)
        thread = allthreads[0]
        thread.current_type = ['prod']
        thread.current_name = ['prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        job_index = 0
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]
        this_inpcrd = ['../test_data/test_velocities.rst7']
        inp_kwargs = {'name': thread.name + '_' + name,
                         'nodes': eval('settings.' + type + '_nodes'),
                         'taskspernode': eval('settings.' + type + '_ppn'),
                         'walltime': eval('settings.' + type + '_walltime'),
                         'mem': eval('settings.' + type + '_mem'),
                         'solver': eval('settings.' + type + '_solver'),
                         'out': thread.name + '_' + name + '.out',
                         'prmtop': thread.topology,
                         'inpcrd': this_inpcrd[job_index],
                         'rst': thread.name + '_' + name + '.rst7',
                         'nc': thread.name + '_' + name + '.nc',
                         'working_directory': settings.working_directory,
                         'extra': eval('settings.' + type + '_extra')}
        numpy.random.seed(0)
        inp = jobtype.get_input_file(thread, job_index, settings, **inp_kwargs)
        assert filecmp.cmp(inp, '../test_data/cp2k_comana_test_velocities.rst7_0_prod.inp')

    def test_get_input_file_cp2k_umbrella_sampling(self):
        """Tests get_input_file with md_engine = 'cp2k' and job_type = 'umbrella_sampling'"""
        settings = configure('../../data/atesa.config')
        settings.md_engine = 'cp2k'
        settings.job_type = 'umbrella_sampling'
        settings.degeneracy = 1
        settings.initial_coordinates = ['../test_data/test_velocities.rst7']
        settings.cvs = ['pytraj.distance(traj, \'@1 @2\')[0]', '(mdtraj.compute_distances(mtraj, numpy.array([[7177, 7178]]))[0][0] * 10) - (mdtraj.compute_distances(mtraj, numpy.array([[4272, 7178]]))[0][0] * 10)']
        settings.us_rc_min = -2
        settings.us_rc_max = 2
        settings.us_rc_step = 0.5
        settings.us_degeneracy = 2
        settings.us_implementation = 'plumed'
        settings.job_type = 'umbrella_sampling'
        settings.as_out_file = '../test_data/as.out'
        allthreads = main.init_threads(settings)
        thread = allthreads[0]
        thread.current_type = ['prod']
        thread.current_name = ['prod']
        jobtype = factory.jobtype_factory(settings.job_type)
        job_index = 0
        type = thread.current_type[job_index]
        name = thread.current_name[job_index]
        this_inpcrd = ['../test_data/test_velocities.rst7']
        thread.history.window = 2.5
        thread.history.index = 1
        inp_kwargs = {'name': thread.name + '_' + name,
                         'nodes': eval('settings.' + type + '_nodes'),
                         'taskspernode': eval('settings.' + type + '_ppn'),
                         'walltime': eval('settings.' + type + '_walltime'),
                         'mem': eval('settings.' + type + '_mem'),
                         'solver': eval('settings.' + type + '_solver'),
                         'out': thread.name + '_' + name + '.out',
                         'prmtop': thread.topology,
                         'inpcrd': this_inpcrd[job_index],
                         'rst': thread.name + '_' + name + '.rst7',
                         'nc': thread.name + '_' + name + '.nc',
                         'working_directory': settings.working_directory,
                         'extra': eval('settings.' + type + '_extra')}
        numpy.random.seed(0)
        inp = jobtype.get_input_file(thread, job_index, settings, **inp_kwargs)
        assert filecmp.cmp(inp, '../test_data/cp2k_umbrella_sampling_2.5_1.in')

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
