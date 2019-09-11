"""
Interface for JobType objects. New JobTypes can be implemented by constructing a new class that inherits from JobType
and implements its abstract methods.
"""

import abc
import os
import sys
import subprocess
import random
try:
    import utilities
except ModuleNotFoundError:
    import atesa_v2.utilities as utilities

class JobType(abc.ABC):
    """
    Abstract base class for job types.

    Implements methods for all of the job type-specific tasks that ATESA might need.

    """

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

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                frame_to_check = self.get_frame(self.traj_files[job_index], -1, settings)
                if utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it
                    os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == C for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
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
        if self.current_type == ['prod', 'prod']:  # aimless shooting only checks termination after prod steps
            thread_terminate = ''       # todo: are there termination criteria to implement for aimless shooting threads?
            global_terminate = False    # todo: implement information error checking

            if thread_terminate:
                self.status = 'terminated after step ' + str(thread.suffix) + ' due to: ' + thread_terminate
                self.terminated = True
            else:
                self.status = 'running step ' + str(thread.suffix + 1)  # suffix isn't updated until call to algorithm()
                self.terminated = False

            return global_terminate

    def update_results(self, allthreads, settings):
        if self.current_type == ['prod', 'prod']:   # aimless shooting only writes after prod steps
            # Initialize as.out if not already extant
            if not os.path.exists(settings.working_directory + '/as.out'):
                open(settings.working_directory + '/as.out', 'w').close()

            # Write CVs to as.out
            open(settings.working_directory + '/as.out', 'a').write(utilities.get_cvs(self.coordinates[0], settings) + '\n')

            # Update current_results, total and accepted move counts, and status.txt
            self.current_results = []   # clear previous results if any
            for job_index in range(len(self.current_type)):
                frame_to_check = self.get_frame(self.traj_files[job_index], -1, settings)
                self.current_results.append(utilities.check_commit(frame_to_check, settings))
                os.remove(frame_to_check)
            self.total_moves += 1
            if self.current_results in [['fwd', 'bwd'], ['bwd', 'fwd']]:
                self.last_accepted = self.traj_files
                self.accept_moves += 1

            with open('status.txt', 'w') as file:
                for thread in allthreads:
                    try:
                        acceptance_percentage = str(100 * thread.accept_moves / thread.total_moves)[0:5] + '%'
                    except ZeroDivisionError:   # 0/0
                        acceptance_percentage = '0%'
                    file.write(thread.initial_coord + ' acceptance ratio: ' + str(thread.accept_moves) + '/' + str(
                            thread.total_moves) + ', or ' + acceptance_percentage + '\n')
                    file.write('  Status: ' + thread.status)
                file.close()

        # Write updated restart.pkl
        pickle.dump(allthreads, open('restart.pkl', 'wb'), protocol=2)  # todo: implement using this file

    def algorithm(self, allthreads, settings):
        # In aimless shooting, algorithm should decide whether or not a new shooting point is needed, obtain it if so,
        # and update self.coordinates to reflect it. Also updates suffix and name attributes.
        if self.current_type == ['prod', 'prod']:
            self.suffix += 1
            self.name = self.initial_coord + '_' + str(self.suffix)
            if self.current_results in [['fwd', 'bwd'], ['bwd', 'fwd']]:    # accepted move
                if random.random() < 0.5:       # randomly select a trajectory (there are only ever two in aimless shooting)
                    job_index = 0
                else:
                    job_index = 1
                frame = random.randint(settings.min_dt, settings.max_dt)
                new_point = self.get_frame(self, self.traj_files[job_index], frame, settings)
                self.coordinates = [new_point]
            else:   # not an accepted move
                if settings.always_new and self.last_accepted:  # choose a new shooting move from last accepted trajectory
                    if random.random() < 0.5:  # randomly select a trajectory (there are only ever two in aimless shooting)
                        job_index = 0
                    else:
                        job_index = 1
                    frame = random.randint(settings.min_dt, settings.max_dt)
                    new_point = self.get_frame(self, self.last_accepted[job_index], frame, settings)
                    self.coordinates = [new_point]
                else:   # always_new = False or there have been no accepted moves in this thread yet
                    pass    # next move will begin from the same point as the last one
        elif self.current_type == ['init']:
            if not os.path.exists(thread.name + '_init.rst7'):  # init step failed, so retry it
                self.current_type = []  # reset current_type so it will be pushed back to ['init'] by thread.process
            else:   # init step produced the desired file, so update coordinates
                self.coordinates = [thread.name + '_init.rst7', utilities.rev_vels(thread.name + '_init.rst7')]


class CommittorAnalysis(JobType):
    """
    Adapter class for committor analysis
    """

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                frame_to_check = self.get_frame(self.traj_files[job_index], -1, settings)
                if utilities.check_commit(frame_to_check, settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it
                    os.remove(frame_to_check)

        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == C for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
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

    def gatekeeper(self, settings):
        # if every job in this thread has status 'C'ompleted/'C'anceled...
        if all(item == C for item in [self.get_status(job_index, settings) for job_index in range(len(self.jobids))]):
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
