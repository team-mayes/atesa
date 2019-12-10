"""
Interface for JobType objects. New JobTypes can be implemented by constructing a new class that inherits from JobType
and implements its abstract methods.
"""

import abc
import os
import sys
import subprocess
import random
import pickle
import argparse
import numpy
import shutil
import time
import pytraj
import warnings
import copy
import re
import psutil
from atesa_v2 import utilities
from atesa_v2 import main
from atesa_v2 import factory

class JobType(abc.ABC):
    """
    Abstract base class for job types.

    Implements methods for all of the job type-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_input_file(self, job_index, settings):
        """
        Obtain appropriate input file for next job.

        At its most simple, implementations of this method can simply return settings.path_to_input_files + '/' +
        settings.job_type + '_' + self.current_type[job_index] + '_' + settings.md_engine + '.in'

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        job_index : int
            0-indexed integer identifying which job within self.current_type to return the input file for
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        input_file : str
            Name of the applicable input file

        """

        pass

    @abc.abstractmethod
    def get_initial_coordinates(self, settings):
        """
        Obtain list of initial coordinate files and copy them to the working directory.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        initial_coordinates : list
            List of strings naming the applicable initial coordinate files that were copied to the working directory

        """

        pass

    @abc.abstractmethod
    def check_for_successful_step(self):
        """
        Check whether a just-completed step was successful, as defined by whether the update_results and
        check_termination methods should be run.

        This method returns True if the previous step appeared successful (distinct from 'accepted' as the term refers
        to for example aimless shooting) and False otherwise. The implementation of what to DO with an unsuccessful
        step should appear in the corresponding algorithm method, which should be run regardless of the output from this
        method.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects

        Returns
        -------
        result : bool
            True if successful; False otherwise

        """

        pass

    @abc.abstractmethod
    def update_history(self, settings, **kwargs):
        """
        Update or initialize the history namespace for this job type.

        This namespace is used to store the full history of a threads coordinate and trajectory files, as well as their
        results if necessary.

        If update_history is called with a kwargs containing {'initialize': True}, it simply prepares a blank
        history namespace and returns it. Otherwise, it adds the values of the desired keywords (which are desired
        depends on the implementation) to the corresponding history attributes in the index given by thread.suffix.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        settings : argparse.Namespace
            Settings namespace object
        kwargs : dict
            Dictionary of arguments that might be used to update the history object

        Returns
        -------
        None

        """

        pass

    @abc.abstractmethod
    def get_inpcrd(self):
        """
        Return a list (possibly of length one) containing the names of the appropriate inpcrd files for the next step
        in the thread given by self.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects

        Returns
        -------
        inpcrd : list
            List of strings containing desired file names

        """

        pass

    @abc.abstractmethod
    def gatekeeper(self, settings):
        """
        Return boolean indicating whether job is ready for next interpretation step.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        status : bool
            If True, ready for next interpretation step; otherwise, False

        """

        pass

    @abc.abstractmethod
    def get_next_step(self, settings):
        """
        Return name of next type of simulation to run in the itinerary of this job type.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        type : list
            List containing strings for name(s) of next type(s) of simulation(s)

        """

        pass

    @abc.abstractmethod
    def get_batch_template(self, type, settings):
        """
        Return name of batch template file for the type of job indicated.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        type : str
            Name of the type of job desired, corresponding to the template file
        settings : argparse.Namespace
                Settings namespace object

        Returns
        -------
        name : str
            Name of the batch file template requested

        """

        pass

    @abc.abstractmethod
    def check_termination(self, allthreads, settings):
        """
        Check termination criteria for the particular thread at hand as well as for the entire process.

        These methods should update self.status and self.terminated in order to communicate the results of the check for
        the inidividual thread, and should return a boolean to communicate the results of the check for the entire run.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        allthreads : list
            The list of all extant Thread objects
        settings : argparse.Namespace
                Settings namespace object

        Returns
        -------
        termination : bool
            Global termination criterion for entire process (True means terminate)

        """

        pass

    @abc.abstractmethod
    def update_results(self, allthreads, settings):
        """
        Update appropriate results file(s) as needed.

        These methods are designed to be called at the end of every step in a job, even if writing to an output file is
        not necessary after that step; for that reason, they also encode the logic to decide when writing is needed.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        allthreads : list
            The list of all extant Thread objects
        settings : argparse.Namespace
                Settings namespace object

        Returns
        -------
        None

        """

        pass

    @abc.abstractmethod
    def algorithm(self, allthreads, running, settings):
        """
        Update thread attributes to prepare for next move.

        This is where the core logical algorithm of a method is implemented. For example, in aimless shooting, it
        encodes the logic of when and how to select a new shooting move.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        allthreads : list
            The list of all extant Thread objects
        running : list
            The list of all currently running Thread objects
        settings : argparse.Namespace
                Settings namespace object

        Returns
        -------
        running : list
            The list of all currently running Thread objects

        """

        pass

    @abc.abstractmethod
    def cleanup(self, settings):
        """
        Perform any last tasks between the end of the main loop and the program exiting.

        Parameters
        ----------
        self : None
            self is not used in cleanup methods
        settings : argparse.Namespace
                Settings namespace object

        Returns
        -------
        None

        """

        pass


class AimlessShooting(JobType):
    """
    Adapter class for aimless shooting
    """

    def get_input_file(self, job_index, settings):
        return settings.path_to_input_files + '/' + settings.job_type + '_' + self.current_type[job_index] + '_' + settings.md_engine + '.in'

    def get_initial_coordinates(self, settings):
        list_to_return = []
        for item in settings.initial_coordinates:
            if settings.degeneracy > 1:     # implements degeneracy option
                og_item = item
                if '/' in item:
                    item = item[item.rindex('/') + 1:]
                list_to_return += [item + '_' + str(this_index) for this_index in range(settings.degeneracy)]
                for file_to_make in list_to_return:
                    shutil.copy(og_item, settings.working_directory + '/' + file_to_make)
            else:
                og_item = item
                if '/' in item:
                    item = item[item.rindex('/') + 1:]
                list_to_return += [item]
                try:
                    shutil.copy(og_item, settings.working_directory + '/' + item)
                except shutil.SameFileError:
                    pass
        return list_to_return

    def check_for_successful_step(self):
        if self.current_type == ['init']:   # requires that self.history.init_coords[-1] exists
            if os.path.exists(self.history.init_coords[-1][0]):
                return True
        if self.current_type == ['prod', 'prod']:   # requires that both files in self.history.prod_trajs[-1] exist
            if all([os.path.exists(self.history.prod_trajs[-1][i]) for i in range(2)]):
                return True
        return False

    def update_history(self, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']:
                self.history = argparse.Namespace()
                self.history.init_inpcrd = []       # list of strings, inpcrd for init steps; initialized by main.init_threads and updated by algorithm
                self.history.init_coords = []       # list of 2-length lists of strings, init [_fwd.rst7, _bwd.rst7]; updated by update_history and then in algorithm
                self.history.prod_trajs = []        # list of 2-length lists of strings, [_fwd.nc, _bwd.nc]; updated by update_history
                self.history.prod_results = []      # list of 2-length lists of strings ['fwd'/'bwd'/'', 'fwd'/'bwd'/'']; updated by update_results
                self.history.last_accepted = -1     # int, index of last accepted prod_trajs entry; updated by update_results (-1 means none yet accepted)
                self.history.timestamps = []        # list of ints representing seconds since the epoch for the end of each step; updated by update_results
            if 'inpcrd' in kwargs.keys():
                self.history.init_inpcrd.append(kwargs['inpcrd'])
        else:   # self.history should already exist
            if self.current_type == ['init']:     # update init attributes
                if 'rst' in kwargs.keys():
                    if len(self.history.init_coords) < self.suffix + 1:
                        self.history.init_coords.append([])
                        if len(self.history.init_coords) < self.suffix + 1:
                            raise IndexError('history.init_coords is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix))
                    self.history.init_coords[self.suffix].append(kwargs['rst'])
            elif self.current_type == ['prod', 'prod']:
                if 'nc' in kwargs.keys():
                    if len(self.history.prod_trajs) < self.suffix + 1:
                        self.history.prod_trajs.append([])
                        if len(self.history.prod_trajs) < self.suffix + 1:
                            raise IndexError('history.prod_trajs is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix))
                    self.history.prod_trajs[self.suffix].append(kwargs['nc'])

    def get_inpcrd(self):
        if self.current_type == ['init']:
            return [self.history.init_inpcrd[-1]]   # should return a list, but history.init_inpcrd contains strings
        elif self.current_type == ['prod', 'prod']:
            return self.history.init_coords[-1]     # should return a list, and history.init_coords contains lists
        else:
            raise ValueError('invalid thread.current_type value: ' + str(self.current_type) + ' for thread: ' + self.history.init_inpcrd[0])

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        if self.current_type == ['prod', 'prod']:
            for job_index in range(len(self.jobids)):
                if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                    frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
                    if frame_to_check and utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
                        self.cancel_job(job_index, settings)        # cancel it
                        os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == 'C' for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
            return True
        else:
            return False

    def get_next_step(self, settings):
        if self.current_type == []:
            self.current_type = ['init']
            self.current_name = ['init']
        elif self.current_type == ['init']:
            self.current_type = ['prod', 'prod']
            self.current_name = ['fwd', 'bwd']
        elif self.current_type == ['prod', 'prod']:
            self.current_type = ['init']
            self.current_name = ['init']
        return self.current_type, self.current_name

    def get_batch_template(self, type, settings):
        if type in ['init', 'prod']:
            templ = settings.md_engine + '_' + settings.batch_system + '.tpl'
            if os.path.exists(settings.path_to_templates + '/' + templ):
                return templ
            else:
                raise FileNotFoundError('cannot find required template file: ' + templ)
        else:
            raise ValueError('unexpected batch template type for aimless_shooting: ' + str(type))

    def check_termination(self, allthreads, settings):
        global_terminate = False    # initialize
        if self.current_type == ['prod', 'prod']:  # aimless shooting only checks termination after prod steps
            thread_terminate = ''
            global_terminate = False

            # Implement settings.max_moves thread termination criterion
            if self.total_moves >= settings.max_moves and settings.max_moves > 0:
                thread_terminate = 'maximum move limit (' + str(settings.max_moves) + ') reached'

            # Implement information error global termination criterion
            # proc_status variable is used to prevent multiple calls to information_error from occuring at once
            if settings.information_error_checking:
                len_data = len(open('as_raw.out', 'r').readlines())     # number of shooting points to date
                if settings.pid == -1:
                    proc_status = 'not_running'     # default value for new process, no process running
                else:
                    try:
                        proc = psutil.Process(settings.pid).status()
                        if proc in [psutil.STATUS_RUNNING, psutil.STATUS_SLEEPING]:
                            proc_status = 'running'
                        elif proc in [psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD]:
                            proc_status = 'not_running'
                        else:
                            warnings.warn('unexpected process state for information_error.py subprocess: ' + proc +
                                          '\nSkipping information error checking at this step')
                            proc_status = 'error'
                    except (psutil.NoSuchProcess, ProcessLookupError):
                        proc_status = 'not_running'
                if (len_data % settings.information_error_freq == 0 and len_data > 0 and proc_status == 'not_running') or (settings.information_error_overdue and proc_status == 'not_running'):
                    # Start separate process calling resample and then information_error
                    process = subprocess.Popen(['resample_and_inferr.py'], stdout=sys.stdout, preexec_fn=os.setsid)
                    settings.pid = process.pid  # store process ID in settings
                    settings.information_error_overdue = False
                elif len_data % settings.information_error_freq == 0 and len_data > 0 and not proc_status == 'not_running':
                    settings.information_error_overdue = True   # so that information error will be checked next time proc_status == 'not_running' regardless of len_data

                # Evaluate termination criterion using KPSS statistic data in info_err.out
                if os.path.exists('info_err.out'):
                    kpss_stat = open('info_err.out', 'r').readlines()[-1].split(' ')[2]     # kpss statistic read from last line
                    if float(kpss_stat) <= 0.05:
                        global_terminate = True     # todo: test information error termination criterion

            if thread_terminate:
                self.status = 'terminated after step ' + str(self.suffix) + ' due to: ' + thread_terminate
                self.terminated = True
            else:
                self.status = 'running step ' + str(self.suffix + 1)  # suffix isn't updated until call to algorithm()
                self.terminated = False

        return global_terminate

    def update_results(self, allthreads, settings):
        if self.current_type == ['prod', 'prod']:   # aimless shooting only writes after prod steps
            # Initialize as.out if not already extant
            if not os.path.exists(settings.working_directory + '/as_raw.out'):
                open(settings.working_directory + '/as_raw.out', 'w').close()

            # Update current_results, total and accepted move counts, and status.txt
            self.history.prod_results.append([])
            for job_index in range(len(self.current_type)):
                frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
                self.history.prod_results[-1].append(utilities.check_commit(frame_to_check, settings))
                os.remove(frame_to_check)
            self.total_moves += 1
            self.history.timestamps.append(int(time.time()))
            if self.history.prod_results[-1] in [['fwd', 'bwd'], ['bwd', 'fwd']]:   # update last accepted move
                if settings.cleanup and self.history.last_accepted >= 0:            # delete previous last accepted trajectories
                    for job_index in range(len(self.current_type)):
                        os.remove(self.history.prod_trajs[self.history.last_accepted][job_index])
                self.history.last_accepted = int(len(self.history.prod_trajs) - 1)   # new index of last accepted move
                self.accept_moves += 1

            # Write CVs to as_raw.out
            if self.history.prod_results[-1][0] in ['fwd', 'bwd']:
                if self.history.prod_results[-1][0] == 'fwd':
                    this_basin = 'B'
                else:   # 'bwd'
                    this_basin = 'A'
                open(settings.working_directory + '/as_raw.out', 'a').write(this_basin + ' <- ')
                open(settings.working_directory + '/as_raw.out', 'a').write(utilities.get_cvs(self.history.init_coords[-1][0], settings) + '\n')
                open(settings.working_directory + '/as_raw.out', 'a').close()

            with open('status.txt', 'w') as file:
                for thread in allthreads:
                    try:
                        acceptance_percentage = str(100 * thread.accept_moves / thread.total_moves)[0:5] + '%'
                    except ZeroDivisionError:   # 0/0
                        acceptance_percentage = '0%'
                    file.write(thread.history.init_inpcrd[0] + ' acceptance ratio: ' + str(thread.accept_moves) +
                               '/' + str(thread.total_moves) + ', or ' + acceptance_percentage + '\n')
                    file.write('  Status: ' + thread.status + '\n')
                file.close()

        # Write updated restart.pkl
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)

    def algorithm(self, allthreads, running, settings):
        # In aimless shooting, algorithm should decide whether or not a new shooting point is needed, obtain it if so,
        # and update self.history to reflect it.
        if self.current_type == ['prod', 'prod']:
            if not all([os.path.exists(self.history.prod_trajs[-1][i]) for i in range(len(self.current_type))]):     # prod step failed, so retry it
                self.current_type = ['init']    # will update to ['prod', 'prod'] in next thread.process call
                return running  # escape immediately

            self.suffix += 1
            self.name = self.history.init_inpcrd[0] + '_' + str(self.suffix)
            if self.history.prod_results[-1] in [['fwd', 'bwd'], ['bwd', 'fwd']]:    # accepted move
                job_index = int(numpy.round(random.random()))    # randomly select a trajectory (there are only ever two in aimless shooting)
                frame = random.randint(settings.min_dt, settings.max_dt)
                new_point = self.get_frame(self.history.prod_trajs[-1][job_index], frame, settings)
                self.history.init_inpcrd.append(new_point)
            else:   # not an accepted move
                if settings.always_new and self.history.last_accepted >= 0:  # choose a new shooting move from last accepted trajectory
                    job_index = int(numpy.round(random.random()))  # randomly select a trajectory (there are only ever two in aimless shooting)
                    frame = random.randint(settings.min_dt, settings.max_dt)
                    new_point = self.get_frame(self.history.prod_trajs[self.history.last_accepted][job_index], frame, settings)
                    self.history.init_inpcrd.append(new_point)
                else:   # always_new = False or there have been no accepted moves in this thread yet
                    self.history.init_inpcrd.append(self.history.init_inpcrd[-1])   # begin next move from same point as last move

            if settings.cleanup and self.history.last_accepted < self.suffix - 1:   # delete this trajectory if it is not the new last_accepted
                for job_index in range(len(self.current_type)):
                    os.remove(self.history.prod_trajs[-1][job_index])

        elif self.current_type == ['init']:
            if not os.path.exists(self.history.init_coords[-1][0]):  # init step failed, so retry it
                self.current_type = []  # reset current_type so it will be pushed back to ['init'] by thread.process
            else:   # init step produced the desired file, so update coordinates
                self.history.init_coords[-1].append(utilities.rev_vels(self.history.init_coords[-1][0]))

        return running

    def cleanup(self, settings):
        # utilities.resample(settings)  # to build as_decorr.out I suppose. Kinda silly.
        pass


# noinspection PyAttributeOutsideInit
class CommittorAnalysis(JobType):
    """
    Adapter class for committor analysis
    """

    def get_input_file(self, job_index, settings):
        return settings.path_to_input_files + '/' + settings.job_type + '_' + self.current_type[job_index] + '_' + settings.md_engine + '.in'

    def get_initial_coordinates(self, settings):
        if settings.committor_analysis_use_rc_out:
            if not os.path.exists(settings.path_to_rc_out):
                raise FileNotFoundError('committor_analysis_use_rc_out = True, but cannot find RC output file at '
                                        'specified path: ' + settings.path_to_rc_out)
            eligible = []       # initialize list of eligible shooting points for committor analysis
            eligible_rcs = []   # another list holding corresponding rc values
            lines = open(settings.path_to_rc_out, 'r').readlines()
            open(settings.path_to_rc_out, 'r').close()
            for line in lines:
                splitline = line.split(': ')             # split line into list [shooting point filename, rc value]
                if abs(float(splitline[1])) <= settings.rc_threshold:
                    eligible.append(splitline[0])
                    eligible_rcs.append(splitline[1])
            if len(eligible) == 0:
                raise RuntimeError('attempted committor analysis, but couldn\'t find any shooting points with reaction '
                                   'coordinate values within ' + str(settings.rc_threshold) + ' of 0 in the RC output '
                                   'file: ' + settings.path_to_rc_out)
            path = settings.path_to_rc_out[:settings.path_to_rc_out.rindex('/')]    # path where rc_out is located
            for item in eligible:
                shutil.copy(path + '/' + item, settings.working_directory + '/' + item)
            return eligible
        else:
            try:
                for item in settings.initial_coordinates:
                    og_item = item
                    if '/' in item:
                        item = item[item.rindex('/') + 1:]
                    try:
                        shutil.copy(og_item, settings.working_directory + '/' + item)
                    except shutil.SameFileError:
                        pass
                return settings.initial_coordinates
            except AttributeError:
                raise RuntimeError('committor_analysis_use_rc_out = False, but initial_coordinates was not provided.')

    def check_for_successful_step(self):
        return True     # nothing to check for in committor analysis

    def update_history(self, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']:
                self.history = argparse.Namespace()
                self.history.prod_inpcrd = []    # one-length list of strings; set by main.init_threads
                self.history.prod_trajs = []     # list of strings; updated by update_history
                self.history.prod_results = []   # list of strings, 'fwd'/'bwd'/''; updated by update_results
            if 'inpcrd' in kwargs.keys():
                self.history.prod_inpcrd.append(kwargs['inpcrd'])
        else:   # self.history should already exist
            if 'nc' in kwargs.keys():
                self.history.prod_trajs.append(kwargs['nc'])

    def get_inpcrd(self):
        return [self.history.prod_inpcrd[0] for null in range(len(self.current_type))]

    def gatekeeper(self, settings):
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
                if frame_to_check and utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it
                    os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == 'C' for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
            return True
        else:
            return False

    def get_next_step(self, settings):
        self.current_type = ['prod' for null in range(settings.committor_analysis_n)]
        self.current_name = [str(int(i)) for i in range(settings.committor_analysis_n)]
        return self.current_type, self.current_name

    def get_batch_template(self, type, settings):
        if type == 'prod':
            templ = settings.md_engine + '_' + settings.batch_system + '.tpl'
            if os.path.exists(settings.path_to_templates + '/' + templ):
                return templ
            else:
                raise FileNotFoundError('cannot find required template file: ' + templ)
        else:
            raise ValueError('unexpected batch template type for committor_analysis: ' + str(type))

    def check_termination(self, allthreads, settings):
        self.terminated = True  # committor analysis threads always terminate after one step
        return False            # no global termination criterion exists for committor analysis

    def update_results(self, allthreads, settings):
        # Initialize committor_analysis.out if not already extant
        if not os.path.exists('committor_analysis.out'):
            open('committor_analysis.out', 'w').write('Committed to Forward Basin / Total Committed to Either Basin')
            open('committor_analysis.out', 'w').close()

        # Update current_results
        for job_index in range(len(self.current_type)):
            frame_to_check = self.get_frame(self.history.prod_trajs[job_index], -1, settings)
            self.history.prod_results.append(utilities.check_commit(frame_to_check, settings))
            os.remove(frame_to_check)

        # Write results to committor_analysis.out
        fwds = 0
        bwds = 0
        for result in self.history.prod_results:
            if result == 'fwd':
                fwds += 1
            elif result == 'bwd':
                bwds += 1
        open('committor_analysis.out', 'a').write(str(int(fwds)) + '/' + str(int(fwds + bwds)) + '\n')

    def algorithm(self, allthreads, running, settings):
        return running    # nothing to set because there is no next step

    def cleanup(self, settings):
        pass


class EquilibriumPathSampling(JobType):
    """
    Adapter class for equilibrium path sampling
    """

    def get_input_file(self, job_index, settings):
        if not settings.eps_n_steps % settings.eps_out_freq == 0:
            raise RuntimeError('eps_n_steps must be evenly divisible by eps_out_freq')

        if self.current_type == ['init']:
            return settings.path_to_input_files + '/' + settings.job_type + '_' + self.current_type[job_index] + '_' + settings.md_engine + '.in'
        else:
            if job_index == 0:  # have to roll to determine the number of fwd and bwd steps
                roll = random.randint(0, int(settings.eps_n_steps/settings.eps_out_freq)) * settings.eps_out_freq  # todo
                self.history.prod_lens.append([roll, settings.eps_n_steps - roll])

            input_file_name = 'eps_' + str(self.history.prod_lens[-1][job_index]) + '.in'

            if not os.path.exists(input_file_name):
                template = settings.env.get_template(settings.md_engine + '_eps_in.tpl')
                filled = template.render(nstlim=str(self.history.prod_lens[-1][job_index]), ntwx=str(settings.eps_out_freq))
                with open(input_file_name, 'w') as newfile:
                    newfile.write(filled)
                    newfile.close()

            return input_file_name

    def get_initial_coordinates(self, settings):
        if settings.as_out_file and settings.rc_reduced_cvs:
            shutil.copy(settings.as_out_file, settings.working_directory + '/')
            if '/' in settings.as_out_file:
                settings.as_out_file = settings.working_directory + '/' + settings.as_out_file[settings.as_out_file.rindex('/') + 1:]
                temp_settings = copy.deepcopy(settings)     # initialize temporary copy of settings to modify
                temp_settings.__dict__.pop('env')           # env attribute is not picklable
                pickle.dump(temp_settings, open(settings.working_directory + '/settings.pkl', 'wb'), protocol=2)
        for item in settings.initial_coordinates:
            og_item = item
            if '/' in item:
                item = item[item.rindex('/') + 1:]
            try:
                shutil.copy(og_item, settings.working_directory + '/' + item)
            except shutil.SameFileError:
                pass
        return settings.initial_coordinates

    def check_for_successful_step(self):
        if self.current_type == ['init']:   # requires that self.history.init_coords[-1] exists
            if os.path.exists(self.history.init_coords[-1][0]):
                return True
        if self.current_type == ['prod', 'prod']:   # requires that both files in self.history.prod_trajs[-1] exist
            if all([os.path.exists(self.history.prod_trajs[-1][i]) for i in range(2)]):
                return True
        return False

    def update_history(self, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']:
                self.history = argparse.Namespace()
                self.history.init_inpcrd = []       # list of strings, inpcrd for init steps; initialized by main.init_threads and updated by algorithm
                self.history.init_coords = []       # list of 2-length lists of strings, init [_fwd.rst7, _bwd.rst7]; updated by update_history and then in algorithm
                self.history.prod_trajs = []        # list of 2-length lists of strings, [_fwd.nc, _bwd.nc]; updated by update_history
                self.history.prod_results = []      # list of 2-length lists of strings ['fwd'/'bwd'/'', 'fwd'/'bwd'/'']; updated by update_results
                self.history.prod_lens = []         # list of 2-length lists of ints indicating the lengths of fwd and bwd trajectories at each step; updated by get_input_file
                self.history.last_accepted = -1     # int, index of last accepted prod_trajs entry; updated by update_results (-1 means none yet accepted)
            if 'inpcrd' in kwargs.keys():
                self.history.init_inpcrd.append(kwargs['inpcrd'])
                cvs = utilities.get_cvs(kwargs['inpcrd'], settings, reduce=settings.rc_reduced_cvs).split(' ')
                init_rc = utilities.evaluate_rc(settings.rc_definition, cvs)
                window_index = 0
                for bounds in settings.eps_bounds:
                    if bounds[0] <= init_rc <= bounds[1]:
                        self.history.bounds = bounds
                        if settings.eps_dynamic_seed:
                            settings.eps_empty_windows[window_index] -= 1  # decrement empty window count in this window
                            if settings.eps_empty_windows[window_index] < 0:  # minimum value 0
                                settings.eps_empty_windows[window_index] = 0
                        break
                    window_index += 1
                try:
                    temp = self.history.bounds  # just to make sure this got set
                except AttributeError:
                    raise RuntimeError('new equilibrium path sampling thread initial coordinates ' + kwargs['inpcrd'] +
                                       ' has out-of-bounds reaction coordinate value: ' + str(init_rc))
        else:   # self.history should already exist
            if self.current_type == ['init']:     # update init attributes
                if 'rst' in kwargs.keys():
                    if len(self.history.init_coords) < self.suffix + 1:
                        self.history.init_coords.append([])
                        if len(self.history.init_coords) < self.suffix + 1:
                            raise IndexError('history.init_coords is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix))
                    self.history.init_coords[self.suffix].append(kwargs['rst'])
            elif self.current_type == ['prod', 'prod']:
                if 'nc' in kwargs.keys():
                    if len(self.history.prod_trajs) < self.suffix + 1:
                        self.history.prod_trajs.append([])
                        if len(self.history.prod_trajs) < self.suffix + 1:
                            raise IndexError('history.prod_trajs is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix))
                    self.history.prod_trajs[self.suffix].append(kwargs['nc'])

    def get_inpcrd(self):
        if self.current_type == ['init']:
            return [self.history.init_inpcrd[-1]]   # should return a list, but history.init_inpcrd contains strings
        elif self.current_type == ['prod', 'prod']:
            return self.history.init_coords[-1]     # should return a list, and history.init_coords contains lists
        else:
            raise ValueError('invalid thread.current_type value: ' + str(self.current_type) + ' for thread: ' + self.history.init_inpcrd[0])

    def gatekeeper(self, settings):
        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == 'C' for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
            return True
        else:
            return False

    def get_next_step(self, settings):
        if self.current_type == []:
            self.current_type = ['init']
            self.current_name = ['init']
        elif self.current_type == ['init']:
            self.current_type = ['prod', 'prod']
            self.current_name = ['fwd', 'bwd']
        elif self.current_type == ['prod', 'prod']:
            self.current_type = ['init']
            self.current_name = ['init']
        return self.current_type, self.current_name

    def get_batch_template(self, type, settings):
        if type in ['init', 'prod']:
            templ = settings.md_engine + '_' + settings.batch_system + '.tpl'
            if os.path.exists(settings.path_to_templates + '/' + templ):
                return templ
            else:
                raise FileNotFoundError('cannot find required template file: ' + templ)
        else:
            raise ValueError('unexpected batch template type for equilibrium path sampling: ' + str(type))

    def check_termination(self, allthreads, settings):
        global_terminate = False    # initialize
        if self.current_type == ['prod', 'prod']:  # equilibrium path sampling only checks termination after prod steps
            thread_terminate = ''
            global_terminate = False    # no global termination criteria for this jobtype

            if settings.samples_per_window > 0:
                samples = 0
                for thread in allthreads:
                    if thread.history.bounds == self.history.bounds:
                        for job_index in range(len(thread.history.prod_results)):
                            samples += len([1 for i in range(len(thread.history.prod_results[job_index])) if thread.history.bounds[0] <= thread.history.prod_results[job_index][i] <= thread.history.bounds[1]])
                if samples >= settings.samples_per_window:
                    thread_terminate = 'reached desired number of samples in this window'

            if thread_terminate:
                self.status = 'terminated after step ' + str(self.suffix) + ' due to: ' + thread_terminate
                self.terminated = True
            else:
                self.status = 'running step ' + str(self.suffix + 1)  # suffix isn't updated until call to algorithm()
                self.terminated = False

        return global_terminate

    def update_results(self, allthreads, settings):
        if self.current_type == ['prod', 'prod']:   # equilibrium path sampling only writes after prod steps
            # Initialize eps.out if not already extant
            if not os.path.exists(settings.working_directory + '/eps.out'):
                open(settings.working_directory + '/eps.out', 'w').write('Lower RC bound; Upper RC bound; RC value')
                open(settings.working_directory + '/eps.out', 'w').close()

            # Update current_results, total and accepted move counts, and status.txt
            self.history.prod_results.append([])
            # First, add results from init frame
            cvs = utilities.get_cvs(self.history.init_coords[-1][0], settings, reduce=settings.rc_reduced_cvs).split(' ')
            self.history.prod_results[-1].append(utilities.evaluate_rc(settings.rc_definition, cvs))
            # Then, add results from fwd and then bwd trajectories, frame-by-frame
            for job_index in range(len(self.current_type)):
                n_frames = pytraj.iterload(self.history.prod_trajs[-1][job_index], settings.topology).n_frames
                for frame in range(n_frames):
                    frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], frame + 1, settings)
                    cvs = utilities.get_cvs(frame_to_check, settings, reduce=settings.rc_reduced_cvs).split(' ')
                    self.history.prod_results[-1].append(utilities.evaluate_rc(settings.rc_definition, cvs))
                    os.remove(frame_to_check)
            self.total_moves += 1
            if True in [self.history.bounds[0] <= rc_value <= self.history.bounds[1] for rc_value in self.history.prod_results[-1]]:
                self.history.last_accepted = int(len(self.history.prod_trajs) - 1)   # new index of last accepted move
                self.accept_moves += 1

            # Write RC values of accepted frames to eps.out
            for rc_value in self.history.prod_results[-1]:
                if self.history.bounds[0] <= rc_value <= self.history.bounds[1]:
                    open(settings.working_directory + '/eps.out', 'a').write(str(self.history.bounds[0]) + ' ' + str(self.history.bounds[1]) + ' ' + str(rc_value) + '\n')
                    open(settings.working_directory + '/eps.out', 'a').close()

            with open('status.txt', 'w') as file:
                for thread in allthreads:
                    try:
                        acceptance_percentage = str(100 * thread.accept_moves / thread.total_moves)[0:5] + '%'
                    except ZeroDivisionError:   # 0/0
                        acceptance_percentage = '0%'
                    file.write(thread.history.init_inpcrd[0] + ' acceptance ratio: ' + str(thread.accept_moves) +
                               '/' + str(thread.total_moves) + ', or ' + acceptance_percentage + '\n')
                    file.write('  Status: ' + thread.status + '\n')
                file.close()

        # Write updated restart.pkl
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)

    def algorithm(self, allthreads, running, settings):
        # In equilibrium path sampling, algorithm should decide whether or not a new shooting point is needed, obtain it
        # if so, and update self.history to reflect it.
        if self.current_type == ['prod', 'prod']:
            self.suffix += 1
            self.name = self.history.init_inpcrd[0] + '_' + str(self.suffix)

            # Need these values for non-accepted move behavior and also for dynamic seeding
            traj = pytraj.iterload(self.history.prod_trajs[self.history.last_accepted][0], settings.topology)
            n_fwd = traj.n_frames
            traj = pytraj.iterload(self.history.prod_trajs[self.history.last_accepted][1], settings.topology)
            n_bwd = traj.n_frames

            if True in [self.history.bounds[0] <= rc_value <= self.history.bounds[1] for rc_value in self.history.prod_results[-1]]:    # accepted move
                traj = pytraj.iterload(self.history.prod_trajs[-1][0], settings.topology)
                n_fwd = traj.n_frames
                traj = pytraj.iterload(self.history.prod_trajs[-1][1], settings.topology)
                n_bwd = traj.n_frames

                if settings.DEBUG:
                    random_bead = n_fwd + n_bwd     # not actually random, for consistency in testing
                else:
                    random_bead = int(random.randint(0, int(n_fwd) + int(n_bwd)))    # randomly select a "bead" from the paired trajectories

                if 0 < random_bead <= n_bwd:
                    frame = random_bead
                    new_point = self.get_frame(self.history.prod_trajs[-1][1], frame, settings)
                elif n_bwd < random_bead <= n_fwd + n_bwd:
                    frame = random_bead - n_bwd
                    new_point = self.get_frame(self.history.prod_trajs[-1][0], frame, settings)
                else:   # random_bead = 0, chooses the "init" bead
                    new_point = self.history.init_inpcrd[-1]
                self.history.init_inpcrd.append(new_point)

            else:   # not an accepted move
                if self.history.last_accepted >= 0:  # choose a new shooting move from last accepted trajectory
                    if settings.DEBUG:
                        random_bead = n_fwd + n_bwd
                    else:
                        random_bead = int(random.randint(0, int(n_fwd) + int(n_bwd)))  # randomly select a "bead" from the paired trajectories

                    if 0 < random_bead <= n_bwd:
                        frame = random_bead
                        new_point = self.get_frame(self.history.prod_trajs[self.history.last_accepted][1], frame, settings)
                    elif n_bwd < random_bead <= n_fwd + n_bwd:
                        frame = random_bead - n_bwd
                        new_point = self.get_frame(self.history.prod_trajs[self.history.last_accepted][0], frame, settings)
                    else:  # random_bead = 0, chooses the "init" bead
                        new_point = self.history.init_inpcrd[self.history.last_accepted]
                    self.history.init_inpcrd.append(new_point)

                else:   # there have been no accepted moves in this thread yet
                    self.history.init_inpcrd.append(self.history.init_inpcrd[-1])   # begin next move from same point as last move

            # Implement dynamic seeding of EPS windows
            if settings.eps_dynamic_seed:
                for rc_value in self.history.prod_results[-1]:  # results are ordered as: [init, fwd, bwd]
                    if not self.history.bounds[0] <= rc_value <= self.history.bounds[1]:
                        try:
                            window_index = [bounds[0] <= rc_value <= bounds[1] for bounds in settings.eps_bounds].index(True)
                        except ValueError:  # rc_value not in range of eps_bounds, so nothing to do here
                            continue
                        if settings.eps_empty_windows[window_index] > 0:    # time to make a new Thread here
                            # First have to get or make the file for the initial coordinates
                            frame_index = self.history.prod_results[-1].index(rc_value)
                            if frame_index == 0:        # init coordinates
                                coord_file = self.history.init_coords[-1][0]
                            elif frame_index <= n_fwd:  # in fwd trajectory
                                coord_file = self.get_frame(self.history.prod_trajs[-1][0], frame_index, settings)
                                if coord_file == '':
                                    raise FileNotFoundError('Trajectory ' + self.history.prod_trajs[-1][0] + ' was not '
                                                            'found in spite of an RC value having been assigned to it.')
                            else:                       # in bwd trajectory
                                coord_file = self.get_frame(self.history.prod_trajs[-1][1], frame_index - n_fwd, settings)
                                if coord_file == '':
                                    raise FileNotFoundError('Trajectory ' + self.history.prod_trajs[-1][1] + ' was not '
                                                            'found in spite of an RC value having been assigned to it.')

                            # Now make the thread and set its parameters
                            new_thread = main.Thread()
                            EquilibriumPathSampling.update_history(new_thread, settings, **{'initialize': True, 'inpcrd': coord_file})
                            new_thread.topology = settings.topology
                            new_thread.name = coord_file + '_' + str(new_thread.suffix)

                            # Append the new thread to allthreads and running
                            allthreads.append(new_thread)
                            running.append(new_thread)

                            if not settings.DEBUG:  # if DEBUG, keep going to maximize coverage
                                return running  # return after initializing one thread to encourage sampling diversity

        elif self.current_type == ['init']:
            if not os.path.exists(self.history.init_coords[-1][0]):  # init step failed, so retry it
                self.current_type = []  # reset current_type so it will be pushed back to ['init'] by thread.process
            else:   # init step produced the desired file, so update coordinates
                self.history.init_coords[-1].append(utilities.rev_vels(self.history.init_coords[-1][0]))

        return running

    def cleanup(self, settings):
        pass


class FindTS(JobType):
    """
    Adapter class for finding transition state (find TS)
    """

    def get_input_file(self, job_index, settings):
        # In find TS, not only do we want to get the input filename, we also want to store the initial basin and write
        # the restraint file to push it into the other basin. First, get the basin...
        commit = utilities.check_commit(self.history.prod_inpcrd[0], settings)
        if commit == '':
            raise RuntimeError('the coordinates provided during a find_ts run must represent a structure in either the '
                               'fwd or bwd basin, but the coordinate file ' + self.history.prod_inpcrd[0] + ' is '
                               'in neither.\nIf it is a transition state guess, you should set jobtype = '
                               '\'aimless_shooting\' instead to begin aimless shooting.')
        elif commit in ['fwd', 'bwd']:
            self.history.init_basin = commit
        else:
            raise RuntimeError('internal error in utilities.check_commit(); did not return valid output with coordinate'
                               ' file: ' + self.history.prod_inpcrd[0])

        # Then, write the restraint file
        if self.history.init_basin == 'fwd':
            other_basin = settings.commit_bwd
        else:  # == 'bwd', the only other valid option
            other_basin = settings.commit_fwd
        mdengine = factory.mdengine_factory(settings.md_engine)
        input_file = mdengine.write_find_ts_restraint(None, other_basin, settings.path_to_input_files + '/find_ts_' + settings.md_engine + '.in')

        return input_file

    def get_initial_coordinates(self, settings):
        if not len(settings.initial_coordinates) == 1:
            raise RuntimeError('when jobtype = \'find_ts\', initial_coordinates must be of length exactly one')

        for item in settings.initial_coordinates:
            og_item = item
            if '/' in item:
                item = item[item.rindex('/') + 1:]
            try:
                shutil.copy(og_item, settings.working_directory + '/' + item)
            except shutil.SameFileError:
                pass

        return settings.initial_coordinates

    def check_for_successful_step(self):
        return True  # nothing to check for in find TS

    def update_history(self, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']:
                self.history = argparse.Namespace()
                self.history.prod_inpcrd = []   # one-length list of strings; set by main.init_threads
                self.history.prod_trajs = []    # list of strings; updated by update_history
                self.history.prod_result = ''   # string, 'fwd'/'bwd'/''; set by update_results or gatekeeper
                self.history.init_basin = ''    # string, 'fwd'/'bwd'/''; set by get_inpcrd
            if 'inpcrd' in kwargs.keys():
                self.history.prod_inpcrd.append(kwargs['inpcrd'])
        else:  # self.history should already exist
            if 'nc' in kwargs.keys():
                self.history.prod_trajs.append(kwargs['nc'])

    def get_inpcrd(self):
        return [self.history.prod_inpcrd[0]]

    def gatekeeper(self, settings):
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':  # if the job in question is running
                frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
                if self.history.init_basin == 'fwd':
                    commit_basin = 'bwd'
                else:   # self.history.init_basin = 'bwd'; get_inpcrd only allows for these two values
                    commit_basin = 'fwd'
                if frame_to_check and utilities.check_commit(frame_to_check, settings) == commit_basin:  # if it has committed to the opposite basin
                    self.history.prod_result = commit_basin
                    self.cancel_job(job_index, settings)  # cancel it
                    os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == 'C' for item in
               [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
            return True
        else:
            return False

    def get_next_step(self, settings):
        self.current_type = ['prod']
        self.current_name = ['ts_guess']
        return self.current_type, self.current_name

    def get_batch_template(self, type, settings):
        if type == 'prod':
            templ = settings.md_engine + '_' + settings.batch_system + '.tpl'
            if os.path.exists(settings.path_to_templates + '/' + templ):
                return templ
            else:
                raise FileNotFoundError('cannot find required template file: ' + templ)
        else:
            raise ValueError('unexpected batch template type for find_ts: ' + str(type))

    def check_termination(self, allthreads, settings):
        self.terminated = True  # find TS threads always terminate after one step
        return False            # no global termination criterion exists for find TS

    # todo: execute the entire find TS process from the end of the simulation to the output right here
    def update_results(self, allthreads, settings):
        # First, store commitment basin if that didn't get set in gatekeeper
        if not self.history.prod_result:
            frame_to_check = self.get_frame(self.history.prod_trajs[-1], -1, settings)
            self.history.prod_result = utilities.check_commit(frame_to_check, settings)

        if self.history.init_basin == 'fwd':
            dir_to_check = 'bwd'
            other_basin_define = settings.commit_bwd
        else:  # self.history.init_basin = 'bwd'; get_inpcrd only allows for these two values
            dir_to_check = 'fwd'
            other_basin_define = settings.commit_fwd

        if not self.history.prod_result == dir_to_check:
            raise RuntimeError('find TS failed to commit to the opposite basin from its given initial coordinates. The '
                               'most likely explanation is that the definition of the target basin ('
                               + dir_to_check + ') is unsuitable for some reason. Please check the produced '
                               'trajectory and output file in the working directory (' +
                               settings.working_directory + ') and either modify the basin definition(s) or the find_ts'
                               ' input file (' + settings.path_to_input_files + '/find_ts_' + settings.md_engine +
                               '.in) accordingly.')

        # Now harvest TSs by inspecting values of target-basin-defining bond lengths vs. frame number

        # We need two helper functions to optimize for the portion of the trajectory most likely to contain the TS
        # There is no need for these functions to be performance-optimized; they will only be used sparingly

        # First just a numerical integrator
        def my_integral(params, my_list):
            # Evaluate a numerical integral of the data in list my_list on the range (params[0],params[1]), where the
            # values in params should be integers referring to indices of my_list.
            partial_list = my_list[int(params[0]):int(params[1])]
            return numpy.trapz(partial_list)

        # And then an optimizer to find the bounds (params) of my_integral that maximize its output
        def my_bounds_opt(objective, data):
            list_of_bounds = []
            for left_bound in range(len(data) - 1):
                for right_bound in range(left_bound + 1, len(data) + 1):
                    if right_bound - left_bound > 1:
                        list_of_bounds.append([left_bound, right_bound])
            output = argparse.Namespace()
            output.best_max = objective(list_of_bounds[0], data)
            output.best_bounds = list_of_bounds[0]
            for bounds in list_of_bounds:
                this_result = objective(bounds, data)
                if this_result > output.best_max:
                    output.best_max = this_result
                    output.best_bounds = bounds
            return output

        # We'll also want a helper function for writing the potential transition state frames to individual files
        def write_ts_guess_frames(traj, frame_indices):
            ts_guesses = []  # initialize list of transition state guesses to test
            for frame_index in frame_indices:
                pytraj.write_traj(self.name + '_ts_guess_' + str(frame_index) + '.rst7', traj,
                                  frame_indices=[frame_index], options='multi', overwrite=True)
                try:
                    os.rename(self.name + '_ts_guess_' + str(frame_index) + '.rst7.1',
                              self.name + '_ts_guess_' + str(frame_index) + '.rst7')
                except FileNotFoundError:
                    if not os.path.exists(self.name + '_ts_guess_' + str(frame_index) + '.rst7'):
                        raise RuntimeError(
                            'Error: attempted to write file ' + self.name + '_ts_guess_' + str(frame_index)
                            + '.rst7, but was unable to. Please ensure that you have '
                              'permission to write to the working directory: ' + settings.working_directory)
                ts_guesses.append(self.name + '_ts_guess_' + str(frame_index) + '.rst7')

            return ts_guesses

        # Now we measure the relevant bond lengths for each frame in the trajectory
        traj = pytraj.iterload(self.history.prod_trajs[0], settings.topology)
        this_lengths = []
        for def_index in range(len(other_basin_define[0])):
            # Each iteration appends a list of the appropriate distances to this_lengths
            this_lengths.append(pytraj.distance(traj, mask='@' + str(other_basin_define[0][def_index]) +
                                                           ' @' + str(other_basin_define[1][def_index])))

        # Now look for the TS by identifying the region in the trajectory with all of the bond lengths at intermediate
        # values (which we'll define as 0.25 < X < 0.75 on a scale of 0 to 1), preferably for several frames in a row.
        norm_lengths = []       # initialize list of normalized (between 0 and 1) lengths
        scored_lengths = []     # initialize list of scores based on lengths
        for lengths in this_lengths:
            normd = [(this_len - min(lengths)) / (max(lengths) - min(lengths)) for this_len in lengths]
            norm_lengths.append(normd)
            scored_lengths.append([(-(this_len - 0.5) ** 2 + 0.0625) for this_len in normd])  # scored on parabola
        sum_of_scores = []
        for frame_index in range(len(scored_lengths[0])):  # sum together scores from each distance at each frame
            sum_of_scores.append(sum([[x[i] for x in scored_lengths] for i in range(len(scored_lengths[0]))][frame_index]))

        # Now I want to find the boundaries between which the TS is most likely to reside. To do this I'll perform a 2D
        # optimization on the integral to find the continuous region that is most positively scored
        opt_result = my_bounds_opt(my_integral, sum_of_scores)

        # If all of sum_of_scores is negative we still want to find a workable max. Also, if there are fewer than five
        # candidate frames we want to test more than that. Either way, shift scores up by 0.1 and try again:
        while opt_result.best_max <= 0 or int(opt_result.best_bounds[1] - opt_result.best_bounds[0]) < 5:
            sum_of_scores = [(item + 0.1) for item in sum_of_scores]
            opt_result = my_bounds_opt(my_integral, sum_of_scores)

        # opt_result.best_bounds now contains the frame indices bounding the region of interest. I want to use pytraj to
        # extract each of these frames as a .rst7 coordinate file and return their names. I'll put a cap on the number
        # of frames that this is done for at, say, 50, so as to avoid somehow ending up with a huge number of candidate
        # TS structures to test.

        if int(opt_result.best_bounds[1] - opt_result.best_bounds[0] + 1) > 50:
            # The best-fit for evenly spacing 50 frames from the bounded region
            frame_indices = [int(ii) for ii in numpy.linspace(opt_result.best_bounds[0], opt_result.best_bounds[1], 50)]
        else:   # number of frames in bounds <= 50, so we'll take all of them
            frame_indices = [int(ii) for ii in range(opt_result.best_bounds[0], opt_result.best_bounds[1]+1)]

        ts_guesses = write_ts_guess_frames(traj, frame_indices)

        ### Here, we do something pretty weird. We want to test the transition state guesses in the ts_guesses list
        ### using aimless shooting, but this isn't an aimless shooting job. So we'll write a new settings object using
        ### the extant one for a template and call main.main() with settings.jobtype = 'aimless_shooting'. Then, we'll
        ### take back control here and use the resulting files to decide what to do next.

        # First, set up the new settings object
        as_settings = settings
        as_settings.initial_coordinates = ts_guesses
        as_settings.working_directory += '/as_test'     # run in a new subdirectory of the current working directory
        as_settings.jobtype = 'aimless_shooting'
        as_settings.max_moves = 10
        as_settings.restart = False
        as_settings.overwrite = True    # in case this is a repeat attempt
        as_settings.cvs = []            # don't need any CVs for this, but this needs to be set to something
        as_settings.resample = False    # just to be safe
        as_settings.information_error_checking = False

        # Because one of the branches below involves repeating this step an unknown (but finite!) number of times, we'll
        # put the whole thing in a while loop controlled by a simple boolean variable.
        not_finished = True

        while not_finished and as_settings.initial_coordinates:     # second condition to ensure initial_coordinates exists
            # Now call main.main() with the new settings object
            main.main(as_settings)   # todo: maybe clean this up sometime for intelligibility
            os.chdir(settings.working_directory)    # chdir back to original working directory

            # Decide what to do next based on the contents of the status.txt file in the new working directory
            if not os.path.exists(as_settings.working_directory + '/status.txt'):
                raise FileNotFoundError('attempted to test transition state guesses in working directory ' +
                                        as_settings.working_directory + ', but the aimless shooting job ended without '
                                        'producing a status.txt file')
            status_lines = open(as_settings.working_directory + '/status.txt', 'r').readlines()
            pattern = re.compile('[0-9.]+\%')
            ratios = []
            names = []
            for line in status_lines:
                if '%' in line:
                    accept_ratio = float(pattern.findall(line)[0].replace('%', ''))/100     # acceptance as a fraction
                    thread_name = line[0:line.index(' ')]
                    ratios.append(accept_ratio)
                    names.append(thread_name)

            if True in [ratio > 0 for ratio in ratios]:     # successful run, so return results and exit
                final_names = []
                for line_index in range(len(ratios)):
                    if ratios[line_index] > 0:
                        final_names.append(names[line_index])
                print('Found working transition state guess(es):\n' + '\n'.join(final_names) + '\nSee ' +
                      as_settings.working_directory + '/status.txt for more information')
                not_finished = False
                continue
            else:       # not so straightforward; this case requires diagnosis
                # Our primary diagnostic tool: the thread.history attributes from the aimless shooting threads
                as_threads = pickle.load(open(as_settings.working_directory + '/restart.pkl', 'rb'))

                # First option: examples of fwd and bwd results found within the same thread indicate that it is a TS,
                # if a poor one; return it and exit
                final_names = []
                for thread in as_threads:
                    reformatted_results = ' '.join([res[0] + ' ' + res[1] for res in thread.history.prod_results])
                    if 'fwd' in reformatted_results and 'bwd' in reformatted_results:
                        final_names.append(thread.history.init_inpcrd[0])   # to match format in status.txt and be maximally useful
                if final_names:
                    print('Found working transition state guess(es):\n' + '\n'.join(final_names) + '\nSee ' +
                          as_settings.working_directory + '/status.txt for more information (acceptance probabilities '
                          'very low; consider trying again with revised basin definitions for better results)')
                    not_finished = False
                    continue

                # Next option: there are only examples of commitment in a single direction. In this case, we want to
                # move the window of test frames over in the appropriate direction and run aimless shooting again
                reformatted_results = ''
                for thread in as_threads:
                    reformatted_results += ' '.join([res[0] + ' ' + res[1] for res in thread.history.prod_results]) + ' '

                if 'fwd' in reformatted_results and not 'bwd' in reformatted_results:       # went to 'fwd' only
                    ts_guesses = []  # re-initialize list of transition state guesses to test
                    if dir_to_check == 'fwd':   # restraint target was 'fwd' so shift to the left (earlier frames)
                        frame_indices = [frame_index - (max(frame_indices) - min(frame_indices)) - 1 for frame_index in frame_indices if frame_index - (max(frame_indices) - min(frame_indices)) - 1 >= 0]
                    if dir_to_check == 'bwd':   # restraint target was 'bwd' so shift to the right (later frames)
                        frame_indices = [frame_index + (max(frame_indices) - min(frame_indices)) + 1 for frame_index in frame_indices if frame_index + (max(frame_indices) - min(frame_indices)) + 1 < traj.n_frames]

                    ts_guesses = write_ts_guess_frames(traj, frame_indices)

                elif 'bwd' in reformatted_results and not 'fwd' in reformatted_results:     # went to 'bwd' only
                    ts_guesses = []  # re-initialize list of transition state guesses to test
                    if dir_to_check == 'bwd':   # restraint target was 'bwd' so shift to the left (earlier frames)
                        frame_indices = [frame_index - (max(frame_indices) - min(frame_indices)) - 1 for frame_index in frame_indices if frame_index - (max(frame_indices) - min(frame_indices)) - 1 >= 0]
                    if dir_to_check == 'fwd':   # restraint target was 'fwd' so shift to the right (later frames)
                        frame_indices = [frame_index + (max(frame_indices) - min(frame_indices)) + 1 for frame_index in frame_indices if frame_index + (max(frame_indices) - min(frame_indices)) + 1 < traj.n_frames]

                    ts_guesses = write_ts_guess_frames(traj, frame_indices)

                # Final option: two possibilities remain here...
                else:
                    if not 'fwd' in reformatted_results and not 'bwd' in reformatted_results:   # no commitment at all
                        raise RuntimeError('no commitments to either basin were observed during aimless shooting.\n'
                                           'The restraint in the barrier crossing simulation appears to have broken the'
                                           ' system in some way, so either try changing your target commitment basin '
                                           'definition and restarting, or you’ll have to use another method to obtain '
                                           'an initial transition state.')
                    else:   # 'fwd' and 'bwd' were both found, but not within a single thread
                        raise RuntimeError('commitments to both the forward and backward basins were observed, but not '
                                           'from any single transition state candidate.\nThe most likely explanation is'
                                           ' that the transition state region is very narrow and lies between two '
                                           'frames in the barrier crossing simulation. Try increasing the frame output '
                                           'frequency or reducing the step size in ' + settings.path_to_input_files +
                                           '/find_ts_' + settings.md_engine + '.in. Alternatively, try adjusting the '
                                           'basin definitions in the config file.')

                # This block only reached if there were only commitments in one direction
                as_settings.initial_coordinates = ts_guesses
                continue    # not_finished = True, so loop will repeat with new initial_coordinates

        if not_finished:    # loop exited because no further ts_guesses were available to test
            raise RuntimeError('transition state guessing was unable to identify a working transition state due to even'
                               ' extreme frames in the barrier crossing simulation failing to commit to the desired '
                               'basin during aimless shooting.\nThe most likely explanation is that one or both basin '
                               'definitions are too loose (states that are not yet truly committed to one or both '
                               'basins are being wrongly identified as committed) or are simply wrong.')

    def algorithm(self, allthreads, running, settings):
        return running  # nothing to set because there is no next step

    def cleanup(self, settings):
        pass
