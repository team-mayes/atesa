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
from atesa_v2 import utilities

class JobType(abc.ABC):
    """
    Abstract base class for job types.

    Implements methods for all of the job type-specific tasks that ATESA might need.

    """

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
    def update_history(self, **kwargs):
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
    def algorithm(self, allthreads, settings):
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

    def check_for_successful_step(self):
        if self.current_type == ['init']:   # requires that self.history.init_coords[-1] exists
            if os.path.exists(self.history.init_coords[-1]):
                return True
        if self.current_type == ['prod', 'prod']:   # requires that both files in self.history.prod_trajs[-1] exist
            if all([os.path.exists(self.history.prod_trajs[-1][i]) for i in range(2)]):
                return True
        return False

    def update_history(self, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize'] == True:
                self.history = argparse.Namespace()
                self.history.init_inpcrd = []    # list of strings, inpcrd for init steps; updated by algorithm
                self.history.init_coords = []    # list of 2-length lists of strings, init [_fwd.rst7, _bwd.rst7]; updated by update_history and then in algorithm
                self.history.prod_trajs = []     # list of 2-length lists of strings, [_fwd.nc, _bwd.nc]; updated by update_history
                self.history.prod_results = []   # list of 2-length lists of strings ['fwd'/'bwd'/'', 'fwd'/'bwd'/'']; updated by update_results
                self.history.last_accepted = -1  # int, index of last accepted prod_trajs entry; updated by update_results (-1 means none yet accepted)
        else:   # self.history should already exist
            if self.current_type == ['init']:     # update init attributes
                if 'rst' in kwargs.keys():
                    if len(self.history.init_coords) < self.suffix + 1:
                        self.history.init_coords.append([])
                        if len(self.history.init_coords) < self.suffix + 1:
                            pickle.dump(allthreads, open('debug.pkl', 'wb'), protocol=2)
                            raise IndexError('history.init_coords is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix) + ''
                                             '\ndumped all threads to debug.pkl')
                    self.history.init_coords[self.suffix].append(kwargs['rst'])
            elif self.current_type == ['prod', 'prod']:
                if 'nc' in kwargs.keys():
                    if len(self.history.prod_trajs) < self.suffix + 1:
                        self.history.prod_trajs.append([])
                        if len(self.history.prod_trajs) < self.suffix + 1:
                            pickle.dump(allthreads, open('debug.pkl', 'wb'), protocol=2)
                            raise IndexError('history.prod_trajs is the wrong length for thread: ' + self.history.init_inpcrd[0] +
                                             '\nexpected length ' + str(self.suffix) + ''
                                             '\ndumped all threads to debug.pkl')
                    self.history.prod_trajs[self.suffix].append(kwargs['nc'])

    def get_inpcrd(self):
        if self.current_type == ['init']:
            return self.history.init_inpcrd
        elif self.current_type == ['prod', 'prod']:
            return self.history.init_coords[-1]
        else:
            pickle.dump(allthreads, open('debug.pkl', 'wb'), protocol=2)
            raise ValueError('invalid thread.current_type value: ' + str(self.current_type) + ' for thread: ' + self.history.init_inpcrd[0] +
                                             '\ndumped all threads to debug.pkl')

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                frame_to_check = self.get_frame(self.history.prod_trajs[job_index], -1, settings)
                if utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
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
        if type in ['init', 'fwd', 'bwd']:
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
            thread_terminate = ''       # todo: are there termination criteria to implement for aimless shooting threads?
            global_terminate = False    # todo: implement information error checking

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
            if not os.path.exists(settings.working_directory + '/as.out'):
                open(settings.working_directory + '/as.out', 'w').close()

            # Write CVs to as.out
            open(settings.working_directory + '/as.out', 'a').write(utilities.get_cvs(self.history.init_coords[-1][0], settings) + '\n')

            # Update current_results, total and accepted move counts, and status.txt
            self.history.prod_results.append([])
            for job_index in range(len(self.current_type)):
                frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
                self.history.prod_results[-1].append(utilities.check_commit(frame_to_check, settings))
                os.remove(frame_to_check)
            self.total_moves += 1
            if self.history.prod_results[-1] in [['fwd', 'bwd'], ['bwd', 'fwd']]:
                self.history.last_accepted = int(len(self.history.prod_trajs) - 1)   # new index of last accepted move
                self.accept_moves += 1

            with open('status.txt', 'w') as file:
                for thread in allthreads:
                    try:
                        acceptance_percentage = str(100 * thread.accept_moves / thread.total_moves)[0:5] + '%'
                    except ZeroDivisionError:   # 0/0
                        acceptance_percentage = '0%'
                    file.write(thread.history.init_inpcrd[0] + ' acceptance ratio: ' + str(thread.accept_moves) +
                               '/' + str(thread.total_moves) + ', or ' + acceptance_percentage + '\n')
                    file.write('  Status: ' + thread.status)
                file.close()

        # Write updated restart.pkl
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)  # todo: implement restarting from this file

    def algorithm(self, allthreads, settings):
        # In aimless shooting, algorithm should decide whether or not a new shooting point is needed, obtain it if so,
        # and update self.coordinates to reflect it. Also updates suffix and name attributes.
        if self.current_type == ['prod', 'prod']:
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
        elif self.current_type == ['init']:
            if not os.path.exists(self.history.init_coords[-1][0]):  # init step failed, so retry it
                self.current_type = []  # reset current_type so it will be pushed back to ['init'] by thread.process
            else:   # init step produced the desired file, so update coordinates
                self.history.init_coords[-1].append(utilities.rev_vels(self.history.init_coords[-1][0]))


class CommittorAnalysis(JobType):
    """
    Adapter class for committor analysis
    """

    def check_for_successful_step(self):    # todo: implement
        pass

    def update_history(self, **kwargs):
        pass    # todo: implement

    def get_inpcrd(self):   # todo: implement
        pass

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                frame_to_check = self.get_frame(self.traj_files[job_index], -1, settings)
                if utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it
                    os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == 'C' for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
            return True
        else:
            return False

    def get_next_step(self, settings):
        if self.current_type == []:
            self.current_type = ['prod' for null in range(settings.committor_analysis_n)]
            self.current_name = [str(int(i)) for i in range(settings.committor_analysis_n)]
        else:
            self.current_type = 'terminate'
            self.current_name = []
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

    def check_termination(self, allthreads, settings):  # todo: implement
        pass

    def update_results(self, allthreads, settings):     # todo: implement
        pass

    def algorithm(self, allthreads, settings):          # todo: implement
        pass


class EquilibriumPathSampling(JobType):
    """
    Adapter class for equilibrium path sampling
    """

    def check_for_successful_step(self):    # todo: implement
        pass

    def update_history(self, **kwargs):
        pass    # todo: implement

    def get_inpcrd(self):   # todo: implement
        pass

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
        if type in ['init', 'fwd', 'bwd']:
            templ = settings.md_engine + '_' + settings.batch_system + '.tpl'
            if os.path.exists(settings.path_to_templates + '/' + templ):
                return templ
            else:
                raise FileNotFoundError('cannot find required template file: ' + templ)
        else:
            raise ValueError('unexpected batch template type for equilibrium path sampling: ' + str(type))

    def check_termination(self, allthreads, settings):  # todo: implement
        pass

    def update_results(self, allthreads, settings):     # todo: implement
        pass

    def algorithm(self, allthreads, settings):          # todo: implement
        pass
