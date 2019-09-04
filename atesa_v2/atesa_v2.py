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
        self.coordinates = ''       # filename containing coordinates
        self.topology = ''          # filename containing topology file
        self.jobids = []            # list of jobids associated with the present step of this thread
        self.traj_files = []        # list of trajectory files associated with the present step of this thread
        self.terminated = False     # boolean indicating whether the thread has reached a termination criterion
        self.current_type = []      # list of job types for the present step of this thread

    def process(self, running, settings):
        return process.process(self, running, settings)

    def interpret(self):
        pass    # todo: implement interpret.py

    def gatekeeper(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.gatekeeper(self, settings)

    def get_next_step(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        self.current_type = jobtype.get_next_step(self, settings)
        return self.current_type

    def get_batch_template(self, type, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.get_batch_template(self, type, settings)

    def get_last_frame(self, job_index, settings):
        mdengine = factory.mdengine_factory(settings.md_engine)
        return mdengine.get_last_frame(self, self.traj_files[job_index], settings)

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
        return allthreads

    # If not restart:
    allthreads = []
    for file in settings.initial_coordinates:
        thread = Thread()
        thread.coordinates = file
        thread.topology = settings.topology
        allthreads.append(thread)
    return allthreads


def check_commit(filename, settings):
    """
    Check commitment of coordinate file to either basin.

    Parameters
    ----------
    filename : str
        Name of .rst7-formatted coordinate file to be checked
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    commit_flag : str
        Either 'fwd' or 'bwd' if the coordinates are in the corresponding basin, or '' if in neither

    """

    traj = pytraj.iterload(filename, settings.topology)
    commit_flag = ''    # initialize

    for i in range(len(settings.commit_fwd[2])):
        if settings.commit_fwd[3][i] == 'lt':
            if pytraj.distance(traj, '@' + str(settings.commit_fwd[0][i]) + ' @' + str(settings.commit_fwd[1][i]),
                               n_frames=1)[0] <= settings.commit_fwd[2][i]:
                commit_flag = 'fwd'     # if a committor test is passed, testing moves on to the next one.
            else:
                commit_flag = ''
                break                   # if a committor test is not passed, all testing in this direction fails
        elif settings.commit_fwd[3][i] == 'gt':
            if pytraj.distance(traj, '@' + str(settings.commit_fwd[0][i]) + ' @' + str(settings.commit_fwd[1][i]),
                               n_frames=1)[0] >= settings.commit_fwd[2][i]:
                commit_flag = 'fwd'
            else:
                commit_flag = ''
                break
        else:
            raise ValueError('An incorrect committor definition \"' + settings.commit_fwd[3][i] + '\" was given for '
                             'index ' + str(i) + ' in the \'fwd\' direction.')

    if commit_flag == '':               # only bother checking for bwd commitment if not fwd committed
        for i in range(len(settings.commit_bwd[2])):
            if settings.commit_bwd[3][i] == 'lt':
                if pytraj.distance(traj, '@' + str(settings.commit_bwd[0][i]) + ' @' + str(settings.commit_bwd[1][i]),
                                   n_frames=1)[0] <= settings.commit_bwd[2][i]:
                    commit_flag = 'bwd' # if a committor test is passed, testing moves on to the next one.
                else:
                    commit_flag = ''
                    break               # if a committor test is not passed, all testing in this direction fails
            elif settings.commit_bwd[3][i] == 'gt':
                if pytraj.distance(traj, '@' + str(settings.commit_bwd[0][i]) + ' @' + str(settings.commit_bwd[1][i]),
                                   n_frames=1)[0] >= settings.commit_bwd[2][i]:
                    commit_flag = 'bwd'
                else:
                    commit_flag = ''
                    break
            else:
                raise ValueError('An incorrect committor definition \"' + settings.commit_bwd[3][i] + '\" was given for'
                                 ' index ' + str(i) + ' in the \'bwd\' direction.')

    return commit_flag


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

    # def gate_keeper(thread):
    #     # Stands between calls to thread.process() and thread.interpret() to prevent the latter from being called until
    #     # all the jobs submitted by the former have completed. In the case of aimless shooting and committor analysis,
    #     # it also cancels those jobs once they've committed to a user-defined energetic basin. Returns True if the
    #     # thread is ready for interpretation and false otherwise.
    #     # todo: consider abstracting this (thread.gate_keeper) with implementations for each type of job. The reason not
    #     # todo: to is because this is already short and complete, but maybe I'll want some kind of extensibility later?
    #     # todo: note: I have since abstracted it; this is kept for now in case I decide to go back.
    #
    #     if settings.job_type in ['aimless_shooting', 'committor_analysis']:
    #         for job_index in range(thread.jobids):
    #             if thread.get_status(job_index, settings) == 'R':                   # if the job in question is running
    #                 if check_commit(thread.get_last_frame(job_index, settings), settings):  # if it has committed to a basin
    #                     thread.cancel_job(job_index, settings)                                  # cancel it
    #
    #     # if every job in this thread has status 'C'ompleted/'C'anceled...
    #     if all(item == C for item in [thread.get_status(job_index, settings) for job_index in range(thread.jobids)]):
    #         return True
    #     else:
    #         return False

    termination_criterion = False   # initialize
    running = allthreads            # to be pruned later by thread.process()

    for thread in allthreads:
        running = thread.process(running, settings)
    while (not termination_criterion) and running:
        for thread in running:
            if thread.gatekeeper(settings):
                thread.interpret()
                running = thread.process(running, settings)

    if termination_criterion:
        print('ATESA run exiting normally (global termination criterion met)')
    else:
        print('ATESA run exiting normally (all threads ended individually)')


if __name__ == "__main__":
    # Obtain settings namespace, initialize threads, and move promptly into main.
    settings = configure.configure(sys.argv[1])
    allthreads = init_threads(settings)
    # allthreads[0].process()
    # sys.exit()
    main(allthreads, settings)
