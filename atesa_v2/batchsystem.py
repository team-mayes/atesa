"""
Interface for BatchSystem objects. New BatchSystems can be implemented by constructing a new class that inherits from
BatchSystem and implements its abstract methods.
"""

import abc
import subprocess

class BatchSystem(abc.ABC):
    """
    Abstract base class for HPC cluster batch systems.

    Implements methods for all of the batch system-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_status(self, jobid, settings):  # todo: this should MAYBE be moved to a method of TaskManager
        """
        Query batch system for a status string for the given job.

        Parameters
        ----------
        jobid : str
            The jobid to query
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        status : str
            A one-character status string. Options are: 'R'unning, 'Q'ueued, and 'C'omplete/'C'anceled.

        """

        pass

    @abc.abstractmethod
    def cancel_job(self, jobid, settings):  # todo: this should DEFINITELY be moved to a method of TaskManager
        """
        Cancel the job given by jobid

        Parameters
        ----------
        jobid : str
            The jobid to cancel
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        output : str
            Raw output string from batch system, if any

        """

        pass
    
    @abc.abstractmethod
    def get_submit_command(self):
        """
        Return the appropriate terminal command for submitting a batch job, with '{file}' where the file to submit 
        should be indicated.

        Parameters
        ----------
        None

        Returns
        -------
        output : str
            Appropriate command including '{file}' substring

        """
        
        pass


class AdaptSlurm(BatchSystem):
    """
    Adapter class for Slurm BatchSystem.

    """

    def get_status(self, jobid, settings):
        if settings.DEBUG:
            return 'C'

        command = 'squeue -o %t --job ' + str(jobid)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()     # decode converts from bytes-like to string
        try:
            output = output.split('\n')[1]
        except IndexError:
            output = 'C'        # job isn't in the queue; so it's 'C'omplete
        if output == 'PD':
            output = 'Q'
        elif output == 'CG':
            output = 'C'
        elif output == 'R':
            output = 'R'        # just to be really explicit
        else:
            raise ValueError('Queried Slurm system status for jobid ' + str(jobid) + ' but got unexpected status: ' + output)

        return output

    def cancel_job(self, jobid, settings):
        command = 'scancel ' + str(jobid)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()     # decode converts from bytes-like to string
        return output
    
    def get_submit_command(self):
        return 'sbatch {file}'


class AdaptPBS(BatchSystem):
    """
    Adapter class for PBS/Torque BatchSystem.

    """

    def get_status(self, jobid, settings):
        if settings.DEBUG:
            return 'C'

        command = 'qstat ' + str(jobid)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()     # decode converts from bytes-like to string
        # Some PBS-specific error handling to help handle common issues by simply resubmitting as necessary.
        while 'Pbs Server is currently too busy to service this request. Please retry this request.' in str(output):
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read()
        if 'Bad UID for job execution MSG=user does not exist in server password file' in str(output) or \
           'This stream has already been closed. End of File.' in str(output):
            time.sleep(60)
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read()

        return output   # todo: this isn't fully implemented; needs to return a one-letter code (see AdaptSlurm)

    def cancel_job(self, jobid, settings):
        command = 'qdel ' + str(jobid)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()     # decode converts from bytes-like to string
        while 'Pbs Server is currently too busy to service this request. Please retry this request.' in str(output):
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read()
        return output

    def get_submit_command(self):
        return 'qsub {file}'
