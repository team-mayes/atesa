"""
atesa_v2.py
Version 2 of Aimless Transition Ensemble Sampling and Analysis refactors the code to make it portable, extensible, and 
flexible.

This script handles the primary loop of building and submitting jobs in independent Threads, using the methods thereof 
to execute various interfaced/abstracted commands.
"""

import sys
import pickle
import pytraj
try:
    import factory
except ModuleNotFoundError:
    import atesa_v2.factory as factory
try:
    import process
except ModuleNotFoundError:
    import atesa_v2.process as process
try:
    import configure
except ModuleNotFoundError:
    import atesa_v2.configure as configure


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
        self.coordinates = []       # filenames containing current coordinates
        self.initial_coord = ''     # filename containing first coordinate file for this thread
        self.topology = ''          # filename containing topology file
        self.jobids = []            # list of jobids associated with the present step of this thread
        self.traj_files = []        # list of trajectory files associated with the present step of this thread
        self.terminated = False     # boolean indicating whether the thread has reached a termination criterion
        self.current_type = []      # list of job types for the present step of this thread
        self.current_name = []      # list of job names corresponding to the job types
        self.current_results = []   # results of each job, if applicable
        self.name = ''              # name of current step
        self.suffix = 1             # index of current step
        self.total_moves = 0        # running total of "moves" attributable to this thread
        self.accept_moves = 0       # running total of "accepted" "moves", as defined by JobType.update_results
        self.last_accepted = []     # list of last accepted trajectories (traj_files of last accepted move)
        self.status = ''            # tag for current status of a thread

    def process(self, running, settings):
        return process.process(self, running, settings)

    def interpret(self, allthreads, settings):
        return interpret.interpret(self, allthreads, settings)

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

    if settings.restart == True:
        allthreads = pickle.load('restart.pkl')
        if settings.restart_terminated_threads:
            for thread in allthreads:
                thread.terminated = False
        return allthreads

    # If not restart:
    allthreads = []
    for file in settings.initial_coordinates:
        thread = Thread()
        thread.coordinates = [file]
        thread.initial_coord = file     # same as coordinates for first step
        thread.topology = settings.topology
        thread.status = 'fresh thread'
        thread.name = thread.initial_coord + '_' + str(thread.suffix)
        allthreads.append(thread)
    return allthreads


def main(allthreads, settings):
    """
    Perform the primary loop of building, submitting, monitoring, and analyzing jobs.

    This function works via a loop of calls to thread.process and thread.interpret for each thread that hasn't
    terminated, until either the global termination criterion is met or all the individual threads have completed.

    Parameters
    ----------
    allthreads : list
        List of all Thread objects, including those for which no further tasks are scheduled.
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    termination_criterion = False   # initialize
    running = allthreads            # to be pruned later by thread.process()

    for thread in allthreads:
        running = thread.process(running, settings)
    while (not termination_criterion) and running:
        for thread in running:
            if thread.gatekeeper(settings):
                termination_criterion = thread.interpret(allthreads, settings)
                if termination_criterion:
                    for thread in running:    # todo: should I replace this with something to finish up running jobs and just block submission of new ones?
                        for job_index in range(len(thread.current_type)):
                            thread.cancel_job(job_index, settings)
                        running = []
                    break
                running = thread.process(running, settings)

    if termination_criterion:
        print('ATESA run exiting normally (global termination criterion met)')
    else:
        print('ATESA run exiting normally (all threads ended individually)')


if __name__ == "__main__":
    # Obtain settings namespace, initialize threads, and move promptly into main.
    settings = configure.configure(sys.argv[1])
    os.chdir(settings.working_directory)
    allthreads = init_threads(settings)
    main(allthreads, settings)
