"""Resampling methods for non-aimless shooting jobs (aimless shooting resampling is in utilities.py)"""

import os
import re
import sys
import glob
import copy
import time
import shutil
import pickle
import warnings
import itertools
import numpy as np
from atesa import factory
from atesa import main
from atesa import utilities
from multiprocess import Pool, Manager, get_context

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

    if settings.transmission_coefficient not in [0, 1, 2]:
        raise RuntimeError('transmission_coefficient set to ' + str(settings.transmission_coefficient) + ' but must be '
                           'either 0, 1, or 2.')

    print('DEBUG (timestamp: ' + str(time.time()) + '): entering resample_committor_analysis')
    sys.stdout.flush()

    jobtype = factory.jobtype_factory('committor_analysis')  # set jobtype

    original_dir = os.getcwd()
    os.chdir(settings.working_directory)  # move to working directory containing committor analysis trajectories
    trajs = glob.glob('*.nc')  # assume all .nc trajectories in the working directory are targets

    if settings.transmission_coefficient in [1, 2]:
        try:    # check if this platform has os.sched_getaffinity
            null = os.sched_getaffinity(0)
            avail = True
        except AttributeError:
            avail = False
        if avail and len(os.sched_getaffinity(0)) > 1:
            print('DEBUG (timestamp: ' + str(time.time()) + '): evaluating transmission coefficient with ' +
                  str(len(os.sched_getaffinity(0))) + ' parallel processes')
            sys.stdout.flush()

            # Break down trajs into chunks, place each chunk in its own list in list trajlists
            trajlists = np.array_split(trajs, len(os.sched_getaffinity(0)))

            # Distribute chunks amongst the workers in a Pool
            with get_context("spawn").Pool(len(os.sched_getaffinity(0))) as p:
                p.starmap(compute_transmission_coefficient,
                          zip([trajs for trajs in trajlists],
                              itertools.repeat(settings),
                              [ii for ii in range(len(os.sched_getaffinity(0)))]))
        else:
            compute_transmission_coefficient(trajs, settings)

        print('DEBUG (timestamp: ' + str(time.time()) + '): computing transmission coefficient')
        sys.stdout.flush()

        # Stitch together results, perform analysis and product output files
        # If desired, compute transmission coefficient as a function of step number
        qdot0 = []
        qdot0_hb_of_t = []
        trajs = []
        for partial_pkl in glob.glob('partial_transmission_coefficient_*.pkl'):
            partial_qdot0, partial_qdot0_hb_of_t, partial_trajs = pickle.load(open(partial_pkl, 'rb'))
            qdot0 += list(partial_qdot0)
            qdot0_hb_of_t += list(partial_qdot0_hb_of_t)
            trajs += list(partial_trajs)

        print('DEBUG (timestamp: ' + str(time.time()) + '): qdot0: ' + str(qdot0))
        print('DEBUG (timestamp: ' + str(time.time()) + '): qdot0_hb_of_t: ' + str(qdot0_hb_of_t))
        sys.stdout.flush()

        open('transmission_coefficient.out', 'w').close()
        kappa = []
        for t in range(min([len(lst) for lst in qdot0_hb_of_t])):
            kappa.append(np.mean([lst[t] for lst in qdot0_hb_of_t]) / (0.5 * np.mean([np.abs(item) for item in qdot0])))
            print('DEBUG (timestamp: ' + str(time.time()) + '): current kappa: ' + str(kappa))
            sys.stdout.flush()
        with open('transmission_coefficient.out', 'a') as f:
            f.write('Transmission coefficient computed as: kappa(n) = <qdot0 * heaviside(rc(n))> / (0.5 * <abs(qdot0)>)'
                    ', where < > indicates the average over each committor analysis trajectory.\nSee Peters\' Reaction '
                    'Rate Theory and Rare Events, Elsevier, 1st ed., Ch. 13, pp. 343\n')
            f.write('Transmission coefficient vs. step number (kappa(n)): ' + str(kappa) + '\n')
            f.write('Initial rates of change of reaction coordinate for each trajectory (qdot0): ' + str(qdot0) + '\n')
            f.write(
                'Rates of change mulitplied by heaviside of reaction coordinate for each step (qdot0 * heaviside(rc(n))): ' + str(
                    qdot0_hb_of_t) + '\n')
            f.write('Trajectory names in the same order as the above lists: ' + str(trajs))

    if settings.transmission_coefficient == 2:  # terminate here if 2; otherwise, must be 0 or 1, so continue on
        return ''

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
    if os.path.exists('committor_analysis.out'):
        os.remove('committor_analysis.out')  # delete old committor_analysis.out before writing new one
    for thread in allthreads:
        jobtype.update_results(thread, allthreads, settings)

    os.chdir(original_dir)  # return to directory from which this was called


def compute_transmission_coefficient(trajs, settings, partial_index=0):
    """
    Resample committor analysis results with new settings without using restart.pkl.

    Go to working directory, find every trajectory file, check its commitment (based on current settings), and based on
    its name, combine those results with other trajectories from the same initial coordinates to produce a new
    committor_analysis.out file.

    Parameters
    ----------
    trajs : list
        List of strings containing names of trajectories to consider
    settings : argparse.Namespace
        Settings namespace object
    partial_index : int
        Integer corresponding to chunk of trajs passed in; used for bookkeeping when parallelized

    Returns
    -------
    None

    """

    print('DEBUG (timestamp: ' + str(time.time()) + '): confirmed transmission_coefficient == True')
    sys.stdout.flush()
    # fwd_commits = 0
    # fwd_crosses = 0
    qdot0 = []  # list of initial rates of change of RC for each trajectory
    qdot0_hb_of_t = []  # list of evaluations of qdot0 times heaviside function hb as a function of step number
    actual_trajs = []

    # Skip if it's already been done; this is just for debugging, not for production
    if os.path.exists('partial_transmission_coefficient_' + str(partial_index) + '.pkl'):
        return 0

    for traj in trajs:
        print('DEBUG (timestamp: ' + str(time.time()) + '): evaluating a trajectory...')
        sys.stdout.flush()
        # commit = utilities.check_commit(traj, settings)
        # this_commits = 0
        # if commit == 'fwd':
        #     fwd_commits += 1
        #     this_commits = 1
        temp_settings = copy.copy(settings)
        temp_settings.include_qdot = False  # must be False for frame='all' on next line
        try:
            cvs_list = utilities.get_cvs(traj, temp_settings, reduce=temp_settings.rc_reduced_cvs, frame='all',
                                     only_in_rc=True).split('\n')  # by-frame list of CV values
        except TypeError as e:
            warnings.warn('Encountered an issue when sampling CVs from trajectory ' + traj + ': ' + str(e) + '\nSkipping...')
            continue
        rcs = [utilities.evaluate_rc(temp_settings.rc_definition, cvs.split()) for cvs in cvs_list]  # by-frame list of RC values
        # this_crosses = 0    # to count the number of times the RC crosses from negative to positive
        # for ii in range(len(rcs) - 1):
        #     if rcs[ii] < 0 and rcs[ii + 1] > 0:
        #         this_crosses += 1
        # if rcs[0] > 0:     # to count the initial cross that would've happened to get to the starting point
        #     this_crosses += 1
        # fwd_crosses += this_crosses     # add to the total
        # with open('transmission_coefficient.out', 'a') as f:
        #     f.write(str(traj) + ': ' + str(this_commits) + '/' + str(this_crosses) + '\n')
        #     if traj == trajs[-1]:   # if this is the last trajectory to look at
        #         f.write('Total: ' + str(fwd_commits) + '/' + str(fwd_crosses) + ' = ' + '%.3f' % float(fwd_commits/fwd_crosses))
        qdot0.append(rcs[1] - rcs[0])  # estimate initial rate of change of RC
        qdot0_hb_of_t.append([qdot0[-1] * np.heaviside(rc, 0) for rc in rcs])
        actual_trajs.append(traj)
        print('DEBUG (timestamp: ' + str(time.time()) + '): current qdot0: ' + str(qdot0))
        # print('DEBUG (timestamp: ' + str(time.time()) + '): current qdot0_hb_of_t: ' + str(qdot0_hb_of_t))
        sys.stdout.flush()

    pickle.dump((qdot0, qdot0_hb_of_t, actual_trajs),
                open('partial_transmission_coefficient_' + str(partial_index) + '.pkl', 'wb'))


def resample_umbrella_sampling(settings):
    """
    Resample umbrella sampling results to produce us_full_cvs.out file.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    original_dir = os.getcwd()
    os.chdir(settings.working_directory)  # move to working directory containing committor analysis trajectories
    temp_settings = copy.deepcopy(settings)
    temp_settings.include_qdot = False  # never want to include_qdot in this upcoming call to get_cvs

    trajs = glob.glob('*.nc')  # all .nc trajectories in the working directory are targets
    with open('us_full_cvs.out', 'w') as f:
        for traj in trajs:
            f.write(utilities.get_cvs(traj, temp_settings, False, 'all') + '\n')

    os.chdir(original_dir)  # return to directory from which this was called