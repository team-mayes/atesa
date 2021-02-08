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
import time
import glob
import psutil
import warnings
import itertools
from atesa import configure
from atesa import factory
from atesa import process
from atesa import interpret
from atesa import utilities
from atesa import information_error
from multiprocess import Pool, Manager

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
        return mdengine.get_frame(traj, frame, settings)

    def get_status(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        return batchsystem.get_status(self.jobids[job_index], settings)

    def cancel_job(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        batchsystem.cancel_job(self.jobids[job_index], settings)


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
            if not thread.current_type == []:
                thread.skip_update = True
        if settings.restart_terminated_threads:
            for thread in allthreads:
                thread.terminated = False
        if settings.job_type == 'aimless_shooting' and settings.information_error_checking:
            if os.path.exists(settings.working_directory + '/info_err.out') and len(open(settings.working_directory + '/info_err.out', 'r').readlines()) > 0:
                info_err_lines = open(settings.working_directory + '/info_err.out', 'r').readlines()

                # Resample completely if there's been a change in the number of definitions of CVs, or in the settings
                # for, information_error_max_dims or information_error_lmax_string
                wrong_length = False
                for data_length in [str(line.split()[0]) for line in info_err_lines]:
                    first_line = open(settings.working_directory + '/as_decorr_' + data_length + '.out', 'r').readlines()[0]
                    num_cvs = len(first_line.replace('A <- ', '').replace('B <- ', '').split())
                    if settings.include_qdot:
                        num_cvs = num_cvs / 2
                    if not num_cvs == len(settings.cvs):
                        wrong_length = True
                if (settings.previous_cvs and not settings.previous_cvs == settings.cvs) or \
                        (not settings.previous_information_error_max_dims == settings.information_error_max_dims) or \
                        (not settings.previous_information_error_lmax_string == settings.information_error_lmax_string) or \
                        wrong_length:
                    utilities.resample(settings, partial=False)
                    information_error.main()

                # Resample if info_err.out is improperly formatted (will not run if resample called above)
                if False in [len(info_err_lines[i].split(' ')) == 2 for i in range(0, len(info_err_lines))]:
                    utilities.resample(settings, partial=True)
                    information_error.main()

                # Resample if info_err.out is missing lines (will not run if resample called above)
                len_data = len(open(settings.working_directory + '/as_raw.out', 'r').readlines())
                last_info_err = info_err_lines[-1].split(' ')[0]
                last_breakpoint = len_data - (len_data % settings.information_error_freq)
                if (last_breakpoint > 0 and not int(last_info_err) == int(last_breakpoint)):
                    utilities.resample(settings, partial=True)
                    information_error.main()
        if settings.job_type == 'equilibrium_path_sampling' and settings.eps_dynamic_seed:    # handle dynamic seeding restart behavior
            for thread in allthreads:
                window_index = 0
                for bounds in settings.eps_bounds:
                    if bounds == thread.history.bounds:
                        settings.eps_empty_windows[window_index] -= 1  # decrement empty window count in this window
                        if settings.eps_empty_windows[window_index] < 0:  # minimum value 0
                            settings.eps_empty_windows[window_index] = 0
                        break
                    window_index += 1

        return allthreads

    # If not restart:
    allthreads = []
    jobtype = factory.jobtype_factory(settings.job_type)

    # Set topology properly even if it's given as a path
    og_prmtop = settings.topology
    if '/' in settings.topology:
        settings.topology = settings.topology[settings.topology.rindex('/') + 1:]
    try:
        shutil.copy(og_prmtop, settings.working_directory + '/' + settings.topology)
    except shutil.SameFileError:
        pass

    for file in jobtype.get_initial_coordinates(settings):
        if '/' in file:
            file = file[file.rindex('/') + 1:]          # drop path to file from filename

        thread = Thread()   # initialize the thread object

        thread.topology = settings.topology

        jobtype.update_history(thread, settings, **{'initialize': True, 'inpcrd': file})    # initialize thread.history

        thread.name = file + '_' + str(thread.suffix)
        allthreads.append(thread)

    return allthreads


def handle_loop_exception(attempted_rescue, running, settings):
    """
    Handle attempted rescue of main loop after encountering an exception, or cancellation of jobs if rescue fails.

    Parameters
    ----------
    attempted_rescue : bool
        True if rescue has already been attempted and this function is being called again. Skips attempting rescue again
        and simply cancels all running jobs.
    running : list
        List of Thread objects that are currently running. These are the threads that will be canceled if the ATESA run
        cannot be rescued.
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    # class UnableToRescueException(Exception):
    #     """ Custom exception for closing out ATESA after an unsuccessful rescue attempt """
    #     pass

    # # todo: finish implementing and then uncomment this
    # if not attempted_rescue:
    #     print('Attempting to remove offending thread(s) and rescue the operation...')
    #     verify_outcome_str, verify_outcome_int = verify_threads.main('restart.pkl')
    #     print(verify_outcome_str)
    #
    #     if verify_outcome_int == 1:  # broken thread removed, continue with attempted rescue
    #         # First, set rescue_running equal to running with deleted threads removed
    #         remaining_threads = pickle.load(open('restart.pkl', 'rb'))
    #         rescue_running = [thread for thread in running if thread in remaining_threads]
    #
    #         # If rescue_running == running, removed threads weren't running so rescue fails. Otherwise...
    #         if not rescue_running == running:
    #             # Then, cancel jobs belonging to deleted threads
    #             deleted_threads = list(set(running) - set(rescue_running))
    #             for thread in deleted_threads:
    #                 try:
    #                     for job_index in range(thread.jobids):
    #                         thread.cancel_job(job_index, settings)
    #                 except Exception as little_e:
    #                     print('Encountered exception while attempting to cancel a job: ' + str(little_e) +
    #                           '\nIgnoring and continuing...')
    #
    #             # Finally, resubmit main() with rescue_running list
    #             main(settings, rescue_running=rescue_running)
    #             return None

    # This code reached if return statement above is not
    print('\nCancelling currently running batch jobs belonging to this process in order to '
          'preserve resources.')
    for thread in running:
        try:
            for job_index in range(len(thread.jobids)):
                thread.cancel_job(job_index, settings)
        except Exception as little_e:
            print('\nEncountered exception while attempting to cancel a job: ' + str(little_e) +
                  '\nIgnoring and continuing...')

    raise RuntimeError('Job cancellation complete, ATESA is now shutting down.')

def main(settings, rescue_running=[]):
    """
    Perform the primary loop of building, submitting, monitoring, and analyzing jobs.

    This function works via a loop of calls to thread.process and thread.interpret for each thread that hasn't
    terminated, until either the global termination criterion is met or all the individual threads have completed.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object
    rescue_running : list
        List of threads passed in from handle_loop_exception, containing running threads. If given, setup is skipped and
        the function proceeds directly to the main loop.

    Returns
    -------
    None

    """

    if not rescue_running:
        # Implement resample
        if settings.job_type == 'aimless_shooting' and settings.resample:
            # Store settings object in the working directory for compatibility with analysis/utility scripts
            if not settings.dont_dump:
                temp_settings = copy.deepcopy(settings)  # initialize temporary copy of settings to modify
                temp_settings.__dict__.pop('env')  # env attribute is not picklable
                pickle.dump(temp_settings, open(settings.working_directory + '/settings.pkl', 'wb'), protocol=2)
            # Run resample
            utilities.resample(settings, partial=False, full_cvs=settings.full_cvs)
            if settings.information_error_checking:     # update info_err.out if called for by settings
                information_error.main()
            sys.exit()

        # Make working directory if it does not exist, handling overwrite and restart as needed
        if os.path.exists(settings.working_directory):
            if settings.overwrite and not settings.restart:
                if os.path.exists(settings.working_directory + '/cvs.txt'):     # a kludge to avoid removing cvs.txt
                    if os.path.exists('ATESA_TEMP_CVS.txt'):
                        raise RuntimeError('tried to create temporary file ATESA_TEMP_CVS.txt in directory: ' +
                                           os.getcwd() + ', but it already exists. Please move, delete, or rename it.')
                    shutil.move(settings.working_directory + '/cvs.txt', 'ATESA_TEMP_CVS.txt')
                shutil.rmtree(settings.working_directory)
                os.mkdir(settings.working_directory)
                if os.path.exists('ATESA_TEMP_CVS.txt'):    # continuation of aforementioned kludge
                    shutil.move('ATESA_TEMP_CVS.txt', settings.working_directory + '/cvs.txt')
            elif not settings.restart and glob.glob(settings.working_directory + '/*') == [settings.working_directory + '/cvs.txt']:
                # Occurs when restart = False, overwrite = False, and auto_cvs is used
                pass
            elif not settings.restart:
                raise RuntimeError('Working directory ' + settings.working_directory + ' already exists, but overwrite '
                                   '= False and restart = False. Either change one of these two settings or choose a '
                                   'different working directory.')
        else:
            if not settings.restart:
                os.mkdir(settings.working_directory)
            else:
                raise RuntimeError('Working directory ' + settings.working_directory + ' does not yet exist, but '
                                   'restart = True.')

        # Store settings object in the working directory for compatibility with analysis/utility scripts
        if os.path.exists(settings.working_directory + '/settings.pkl'):    # for checking for need for resample later
            previous_settings = pickle.load(open(settings.working_directory + '/settings.pkl', 'rb'))
            settings.previous_cvs = previous_settings.cvs
            try:
                settings.previous_information_error_max_dims = previous_settings.information_error_max_dims
            except AttributeError:
                pass
            try:
                settings.previous_information_error_lmax_string = previous_settings.information_error_lmax_string
            except AttributeError:
                pass
        if not settings.dont_dump:
            temp_settings = copy.deepcopy(settings)  # initialize temporary copy of settings to modify
            temp_settings.__dict__.pop('env')  # env attribute is not picklable (update: maybe no longer true, but doesn't matter)
            pickle.dump(temp_settings, open(settings.working_directory + '/settings.pkl', 'wb'), protocol=2)

        # Build or load threads
        allthreads = init_threads(settings)

        # Move runtime to working directory
        os.chdir(settings.working_directory)

        running = allthreads.copy()     # to be pruned later by thread.process()
        attempted_rescue = False        # to keep track of general error handling below
    else:
        allthreads = pickle.load(open(settings.working_directory + '/restart.pkl', 'rb'))
        running = rescue_running
        attempted_rescue = True

    # Initialize threads with first process step
    try:
        if not rescue_running:  # if rescue_running, this step has already finished and we just want the while loop
            for thread in allthreads:
                running = thread.process(running, settings)
    except Exception as e:
        if settings.restart:
            print('The following error occurred while attempting to initialize threads from restart.pkl. It may be '
                  'corrupted.')
                  #'If you haven\'t already done so, consider running verify_threads.py to remove corrupted threads from this file.'
        raise e

    try:
        if settings.job_type == 'aimless_shooting' and len(os.sched_getaffinity(0)) > 1:
            # Initialize Manager for shared data across processes; this is necessary because multiprocessing is being
            # retrofitted to code designed for serial processing, but it works!
            manager = Manager()

            # Setup Managed allthreads list
            managed_allthreads = []
            for thread in allthreads:
                thread_dict = thread.__dict__
                thread_history_dict = thread.history.__dict__
                managed_thread = Thread()
                managed_thread.history = manager.Namespace()
                managed_thread.__dict__.update(thread_dict)
                managed_thread.history.__dict__.update(thread_history_dict)
                managed_allthreads.append(managed_thread)
            allthreads = manager.list(managed_allthreads)

            # Setup Managed settings Namespace
            settings_dict = settings.__dict__
            managed_settings = manager.Namespace()
            # Need to explicitly update every key because of how the Managed Namespace works.
            # Calling exec is the best way to do this I could find. Updating managed_settings.__dict__ doesn't work.
            for key in settings_dict.keys():
                exec('managed_settings.' + key + ' = settings_dict[key]')

            # Distribute processes among available core Pool
            with Pool(len(os.sched_getaffinity(0))) as p:
                p.starmap(main_loop, zip(itertools.repeat(managed_settings), itertools.repeat(allthreads), [[thread] for thread in allthreads]))
        else:
            main_loop(settings, allthreads, running)
    except AttributeError:  # os.sched_getaffinity raises AttributeError on non-UNIX systems.
        main_loop(settings, allthreads, running)

    ## Deprecated thread pool
    # pool = ThreadPool(len(allthreads))
    # func = partial(main_loop, settings)
    # results = pool.map(func, [[thread] for thread in allthreads])

    jobtype = factory.jobtype_factory(settings.job_type)
    jobtype.cleanup(settings)

    return 'ATESA run exiting normally'

def main_loop(settings, allthreads, running):
    termination_criterion = False
    attempted_rescue = False
    interpreted = []
    # Begin main loop
    # This whole thing is in a try-except block to handle cancellation of jobs when the code crashes in any way
    try:
        while (not termination_criterion) and running:
            for thread in running:
                if thread.gatekeeper(settings):
                    termination_criterion, running = thread.interpret(allthreads, running, settings)
                    if attempted_rescue == True:
                        interpreted.append(thread)
                    if termination_criterion:
                        for thread in running:
                            for job_index in range(len(thread.current_type)):
                                thread.cancel_job(job_index, settings)
                        running = []
                        if not settings.pid == -1:  # finish up currently running resample_and_inferr, if appropriate
                            proc_status = 'running'
                            while proc_status == 'running':
                                try:
                                    proc = psutil.Process(settings.pid).status()
                                    if proc in [psutil.STATUS_RUNNING, psutil.STATUS_SLEEPING, psutil.STATUS_DISK_SLEEP]:
                                        proc_status = 'running'
                                        time.sleep(60)  # wait 1 minute before checking again
                                    else:
                                        proc_status = 'not_running'
                                except (psutil.NoSuchProcess, ProcessLookupError):
                                    proc_status = 'not_running'
                        break
                    running = thread.process(running, settings)
                else:
                    time.sleep(30)  # to prevent too-frequent calls to batch system by thread.gatekeeper

            if all([thread in interpreted for thread in running]):
                attempted_rescue = False    # every thread has passed at an interpret step, so rescue was successful!

    except Exception as e:
        print(str(e))
        handle_loop_exception(attempted_rescue, running, settings)


def run_main():
    # Obtain settings namespace, initialize threads, and move promptly into main.
    try:
        working_directory = sys.argv[2]
    except IndexError:
        working_directory = ''
    try:
        settings = configure.configure(sys.argv[1], working_directory)
    except IndexError:
        raise RuntimeError('No configuration file specified. See documentation at atesa.readthedocs.io for details.')
    exit_message = main(settings)
    print(exit_message)

if __name__ == "__main__":
    run_main()
