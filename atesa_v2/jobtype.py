"""
Interface for JobType objects. New JobTypes can be implemented by constructing a new class that inherits from JobType
and implements its abstract methods.
"""

import abc
import os
import sys
import subprocess

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


class AimlessShooting(JobType):
    """
    Adapter class for aimless shooting
    """

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                if check_commit(self.get_last_frame(job_index, settings), settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it

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


class CommittorAnalysis(JobType):
    """
    Adapter class for committor analysis
    """

    def gatekeeper(self, settings):
        # Implement flexible length shooting...
        for job_index in range(len(self.jobids)):
            if self.get_status(job_index, settings) == 'R':     # if the job in question is running
                if check_commit(self.get_last_frame(job_index, settings), settings):  # if it has committed to a basin
                    self.cancel_job(job_index, settings)        # cancel it

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
            raise ValueError('unexpected batch template type for equilbrium path sampling: ' + str(type))
