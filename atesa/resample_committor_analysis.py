"""Just a quick routine for resampling committor analysis results with a new settings object and without assuming the
existence of a restart.pkl file in the working directory."""

import os
import re
import glob
import shutil
from atesa import factory
from atesa import main

def resample_committor_analysis(settings):
    """
    Resample committor analysis results with new settings without using restart.pkl.

    Go to working directory, find every trajectory file, check its commitment (based on current settings), and based on
    its name, combine those results with other trajectories from the same initial coordinates to produce a new
    committor_analysis.out file.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """
    original_dir = os.getcwd()
    os.chdir(settings.working_directory)    # move to working directory containing committor analysis trajectories

    jobtype = factory.jobtype_factory('committor_analysis')  # set jobtype

    trajs = glob.glob('*.nc')   # assume all .nc trajectories in the working directory are targets

    # First, make new thread objects based on names of trajectories in trajs. We can take advantage of the fact that all
    # the trajectories are names as: [threadname]_[int]_[another_int].nc
    pattern = re.compile('_[0-9]+_[0-9]+\.nc')
    allthreads = []
    for traj in trajs:
        try:
            threadname = traj.replace(pattern.findall(traj)[-1], '')
        except IndexError:
            raise RuntimeError('Tried to get thread name from trajectory file: ' + traj + ', but it was not formatted '
                               'as expected. Are you sure this is an ATESA committor analysis directory?')

        # Create new thread for this trajectory if necessary
        if not threadname in [thread.name for thread in allthreads]:
            try:
                assert os.path.exists(threadname)   # initial coordinate file for this thread needs to exist
            except AssertionError:
                raise RuntimeError('Found trajectories appearing to belong to a thread based on initial coordinate file'
                                   ' ' + threadname + ' but no such file was found in the working directory. Either '
                                   'this isn\'t an ATESA committor analysis directory or it\'s corrupted or modified.')

            thread = main.Thread()

            # Initialize thread (modified from main.init_threads)
            og_prmtop = settings.topology
            if '/' in settings.topology:
                settings.topology = settings.topology[settings.topology.rindex('/') + 1:]
            try:
                shutil.copy(og_prmtop, settings.working_directory + '/' + settings.topology)
            except shutil.SameFileError:
                pass

            thread.topology = settings.topology
            thread.current_type = []

            jobtype.update_history(thread, settings, **{'initialize': True, 'inpcrd': threadname})  # initialize thread.history

            thread.name = threadname
            allthreads.append(thread)

        # Append this trajectory to the list of moves in this thread
        thread = allthreads[[threadname == thread.name for thread in allthreads].index(True)]
        thread.history.prod_trajs.append(traj)
        thread.current_type.append('prod')

    # Now we can use the extant committor analysis method for getting results for each thread
    for thread in allthreads:
        jobtype.update_results(thread, allthreads, settings)

    os.chdir(original_dir)  # return to directory from which this was called
