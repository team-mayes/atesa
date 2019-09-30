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
import pytraj
from atesa_v2 import utilities

class JobType(abc.ABC):
    """
    Abstract base class for job types.

    Implements methods for all of the job type-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_initial_coordinates(self, settings):
        """
        Obtain list of initial coordinate files.

        At its most simple, implementations of this method can simply return settings.initial_coordinates.

        Parameters
        ----------
        self : Thread
            Methods in the JobType abstract base class are intended to be invoked by Thread objects
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        initial_coordinates : list
            List of strings naming the applicable initial coordinate files

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

    def get_initial_coordinates(self, settings):
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
                self.history.init_inpcrd = []    # list of strings, inpcrd for init steps; initialized by main.init_threads and updated by algorithm
                self.history.init_coords = []    # list of 2-length lists of strings, init [_fwd.rst7, _bwd.rst7]; updated by update_history and then in algorithm
                self.history.prod_trajs = []     # list of 2-length lists of strings, [_fwd.nc, _bwd.nc]; updated by update_history
                self.history.prod_results = []   # list of 2-length lists of strings ['fwd'/'bwd'/'', 'fwd'/'bwd'/'']; updated by update_results
                self.history.last_accepted = -1  # int, index of last accepted prod_trajs entry; updated by update_results (-1 means none yet accepted)
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

            # Write CVs to as.out
            if self.history.prod_results[-1][0] in ['fwd', 'bwd']:
                if self.history.prod_results[-1][0] == 'fwd':
                    this_basin = 'B'
                else:   # 'bwd'
                    this_basin = 'A'
                open(settings.working_directory + '/as.out', 'a').write(this_basin + ' <- ')
                open(settings.working_directory + '/as.out', 'a').write(utilities.get_cvs(self.history.init_coords[-1][0], settings) + '\n')
                open(settings.working_directory + '/as.out', 'a').close()

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

    def algorithm(self, allthreads, settings):
        # In aimless shooting, algorithm should decide whether or not a new shooting point is needed, obtain it if so,
        # and update self.history to reflect it.
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


# noinspection PyAttributeOutsideInit
class CommittorAnalysis(JobType):
    """
    Adapter class for committor analysis
    """

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
            return eligible
        else:
            try:
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

    def algorithm(self, allthreads, settings):
        pass    # nothing to set because there is no next step


class EquilibriumPathSampling(JobType):
    """
    Adapter class for equilibrium path sampling
    """

    def get_initial_coordinates(self, settings):
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
                self.history.init_inpcrd = []    # list of strings, inpcrd for init steps; initialized by main.init_threads and updated by algorithm
                self.history.init_coords = []    # list of 2-length lists of strings, init [_fwd.rst7, _bwd.rst7]; updated by update_history and then in algorithm
                self.history.prod_trajs = []     # list of 2-length lists of strings, [_fwd.nc, _bwd.nc]; updated by update_history
                self.history.prod_results = []   # list of 2-length lists of strings ['fwd'/'bwd'/'', 'fwd'/'bwd'/'']; updated by update_results
                self.history.last_accepted = -1  # int, index of last accepted prod_trajs entry; updated by update_results (-1 means none yet accepted)
            if 'inpcrd' in kwargs.keys():
                self.history.init_inpcrd.append(kwargs['inpcrd'])
                cvs = utilities.get_cvs(kwargs['inpcrd'], settings, reduce=settings.rc_reduced_cvs).split(' ')
                init_rc = utilities.evaluate_rc(settings.rc_definition, cvs)
                for bounds in settings.eps_bounds:
                    if bounds[0] <= init_rc <= bounds[1]:
                        self.history.bounds = bounds
                        break
                try:
                    temp = self.history.bounds  # just to make sure this got set
                except NameError:
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
        if type in ['init', 'fwd', 'bwd']:
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
            thread_terminate = ''       # todo: are there termination criteria to implement for equilibrium path sampling threads?
            global_terminate = False    # no global termination criteria for this jobtype

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
            for job_index in range(len(self.current_type)):
                frame_to_check = self.get_frame(self.history.prod_trajs[-1][job_index], -1, settings)
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

    def algorithm(self, allthreads, settings):
        # In equilibrium path sampling, algorithm should decide whether or not a new shooting point is needed, obtain it
        # if so, and update self.history to reflect it.
        if self.current_type == ['prod', 'prod']:
            self.suffix += 1
            self.name = self.history.init_inpcrd[0] + '_' + str(self.suffix)
            if True in [self.history.bounds[0] <= rc_value <= self.history.bounds[1] for rc_value in self.history.prod_results[-1]]:    # accepted move
                traj = pytraj.iterload(self.history.prod_trajs[-1][0], settings.topology)
                n_fwd = traj.n_frames
                traj = pytraj.iterload(self.history.prod_trajs[-1][1], settings.topology)
                n_bwd = traj.n_frames
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
                    traj = pytraj.iterload(self.history.prod_trajs[self.history.last_accepted][0], settings.topology)
                    n_fwd = traj.n_frames
                    traj = pytraj.iterload(self.history.prod_trajs[self.history.last_accepted][1], settings.topology)
                    n_bwd = traj.n_frames
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
        elif self.current_type == ['init']:
            if not os.path.exists(self.history.init_coords[-1][0]):  # init step failed, so retry it
                self.current_type = []  # reset current_type so it will be pushed back to ['init'] by thread.process
            else:   # init step produced the desired file, so update coordinates
                self.history.init_coords[-1].append(utilities.rev_vels(self.history.init_coords[-1][0]))
