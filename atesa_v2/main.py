"""
main.py
Version 2 of Aimless Transition Ensemble Sampling and Analysis refactors the code to make it portable, extensible, and 
flexible.

This script handles the primary loop of building and submitting jobs in independent Threads, using the methods thereof 
to execute various interfaced/abstracted commands.
"""

import sys
import os
import shutil
import pickle
import pytraj
import copy
from atesa_v2 import configure
from atesa_v2 import factory
from atesa_v2 import process
from atesa_v2 import interpret
from atesa_v2 import utilities

import tracemalloc

class Thread(object):
    """
    Object representing a series of simulations and containing the relevant information to define its current state.

    Threads represent the level on which ATESA is parallelized. This flexible object is used for every type of job
    performed by ATESA.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    def __init__(self):
        self.topology = ''              # filename containing topology file
        self.jobids = []                # list of jobids associated with the present step of this thread
        self.terminated = False         # boolean indicating whether the thread has reached a termination criterion
        self.current_type = []          # list of job types for the present step of this thread
        self.current_name = []          # list of job names corresponding to the job types
        self.current_results = []       # results of each job, if applicable
        self.name = ''                  # name of current step
        self.suffix = 0                 # index of current step
        self.total_moves = 0            # running total of "moves" attributable to this thread
        self.accept_moves = 0           # running total of "accepted" "moves", as defined by JobType.update_results
        self.status = 'fresh thread'    # tag for current status of a thread
        self.skip_update = False        # used by restart to resubmit jobs as they were rather than doing the next step

    def process(self, running, settings):
        return process.process(self, running, settings)

    def interpret(self, allthreads, running, settings):
        return interpret.interpret(self, allthreads, running, settings)

    def gatekeeper(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.gatekeeper(self, settings)

    def get_next_step(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        self.current_type, self.current_name = jobtype.get_next_step(self, settings)
        return self.current_type, self.current_name

    def get_batch_template(self, type, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.get_batch_template(self, type, settings)

    def get_frame(self, traj, frame, settings):
        mdengine = factory.mdengine_factory(settings.md_engine)
        return mdengine.get_frame(self, traj, frame, settings)

    def get_status(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        return batchsystem.get_status(self, self.jobids[job_index], settings)

    def cancel_job(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        batchsystem.cancel_job(self, self.jobids[job_index], settings)


def init_threads(settings):
    """
    Initialize all the Thread objects called for by the user input file.

    In the case where settings.restart == True, this involves unpickling restart.pkl; otherwise, brand new objects are
    produced in accordance with settings.job_type (aimless_shooting, committor_analysis, equilibrium_path_sampling, or
    isee).

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    allthreads : list
        List of all Thread objects, including those for which no further tasks are scheduled.

    """

    if settings.restart:
        allthreads = pickle.load(open(settings.working_directory + '/restart.pkl', 'rb'))
        for thread in allthreads:
            thread.skip_update = True
        if settings.restart_terminated_threads:
            for thread in allthreads:
                thread.terminated = False
        return allthreads

    # If not restart:
    allthreads = []
    jobtype = factory.jobtype_factory(settings.job_type)
    for file in jobtype.get_initial_coordinates(None, settings):
        if '/' in file:
            file = file[file.rindex('/') + 1:]          # drop path to file from filename

        thread = Thread()
        jobtype.update_history(thread, settings, **{'initialize': True, 'inpcrd': file})
        thread.topology = settings.topology
        thread.name = file + '_' + str(thread.suffix)
        allthreads.append(thread)

    return allthreads


def main(settings):
    """
    Perform the primary loop of building, submitting, monitoring, and analyzing jobs.

    This function works via a loop of calls to thread.process and thread.interpret for each thread that hasn't
    terminated, until either the global termination criterion is met or all the individual threads have completed.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    if settings.resample:
        utilities.resample(settings)
        sys.exit()

    # Make working directory if it does not exist, handling overwrite and restart as needed
    if os.path.exists(settings.working_directory):
        if settings.overwrite:
            shutil.rmtree(settings.working_directory)
            os.mkdir(settings.working_directory)
        elif not settings.restart:
            raise RuntimeError('Working directory ' + settings.working_directory + ' already exists, but overwrite = '
                               'False and restart = False. Either change one of these two settings or choose a '
                               'different working directory.')
    else:
        if not settings.restart:
            os.mkdir(settings.working_directory)
        else:
            raise RuntimeError('Working directory ' + settings.working_directory + ' does not yet exist, but restart = '
                               'True.')

    # Build threads and move necessary files to working directory
    allthreads = init_threads(settings)
    shutil.copy(settings.topology, settings.working_directory)

    # Move runtime to working directory
    os.chdir(settings.working_directory)

    # Store settings object in the working directory for posterity and for compatibility with analysis/utility scripts
    temp_settings = copy.deepcopy(settings)        # initialize temporary copy of settings to modify
    temp_settings.__dict__.pop('env')                           # env attribute is not picklable
    pickle.dump(temp_settings, open('settings.pkl', 'wb'), protocol=2)

    termination_criterion = False   # initialize global termination criterion boolean
    running = allthreads            # to be pruned later by thread.process()

    # Begin main loop
    for thread in allthreads:
        running = thread.process(running, settings)
    while (not termination_criterion) and running:

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)

        for thread in running:
            if thread.gatekeeper(settings):
                termination_criterion, running = thread.interpret(allthreads, running, settings)
                if termination_criterion:
                    for thread in running:    # todo: should I replace this with something to finish up running jobs and just block submission of new ones?
                        for job_index in range(len(thread.current_type)):
                            thread.cancel_job(job_index, settings)
                        running = []
                    break
                running = thread.process(running, settings)

    jobtype = factory.jobtype_factory(settings.job_type)
    jobtype.cleanup(None, settings)

    if termination_criterion:
        print('ATESA run exiting normally (global termination criterion met)')
    else:
        print('ATESA run exiting normally (all threads ended individually)')

if __name__ == "__main__":
    tracemalloc.start()
    # Obtain settings namespace, initialize threads, and move promptly into main.
    settings = configure.configure(sys.argv[1]) #'data/atesa.config')
    main(settings)
