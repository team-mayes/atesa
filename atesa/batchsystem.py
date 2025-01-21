"""
Interface for BatchSystem objects. New BatchSystems can be implemented by constructing a new class that inherits from
BatchSystem and implements its abstract methods.
"""

import re
import abc
import subprocess
import psutil
import time

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
    def submit_batch(self, filename, settings):
        """
        Submit a batch file to the task manager.

        Parameters
        ----------
        filename : str
            Name of batch file to submit
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        jobid : str
            Identification number for this task, such that it can be cancelled by referring to this string

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
def get_file_suffix(self):
    """
    Return the appropriate suffix for batch files (not including ".")

    Parameters
    ----------
    None

    Returns
    -------
    output : str
        Appropriate suffix

    """

    pass


class AdaptSlurm(BatchSystem):
    """
    Adapter class for Slurm BatchSystem.

    """

    def get_status(self, jobid, settings):
        if settings.DEBUG:
            return 'C'

        count = 1
        max_tries = 5
        output = 'first_attempt'
        errors = ['first_attempt', 'slurm_load_jobs', 'slurm_receive_msg', 'slurm_msg_sendto', 'send/recv', 'connect failure']   # error messages to retry on

        command = 'squeue -o %t --job ' + str(jobid)
        while True in [error in output for error in errors] and count <= max_tries and not 'Invalid job id specified' in output:
            if not output == 'first_attempt':
                time.sleep(30)      # wait 30 seconds before trying again
                count += 1
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read().decode()     # decode converts from bytes-like to string

        # Raise error if anything other than a status or invalid job id
        if not output.split('\n')[0] == 'ST' and not 'Invalid job id' in output:
            raise ValueError('Queried Slurm system status for jobid ' + str(jobid) + ' but got unexpected status: ' + output)

        try:
            output = output.split('\n')[1]  # get status code; if invalid job id, this is an empty string
        except IndexError:
            output = output     # should never happen but if it does somehow then this makes the error message useful
        if output in ['PD', 'CF']:
            output = 'Q'
        elif output in ['CG', 'ST'] or 'Invalid job id' in output or not output:  # 'Invalid job id' or empty output returned when job is finished in Slurm
            output = 'C'
        elif output == 'R':
            output = 'R'            # just to be really explicit I guess
        else:
            raise ValueError('Queried Slurm system status for jobid ' + str(jobid) + ' but got unexpected status: ' + output)

        return output

    def submit_batch(self, filename, settings):
        command = 'sbatch {file}'.replace('{file}', filename)

        if settings.DEBUG:
            command = ('echo "this is a nonsense string for testing purposes: 123456, '
                       'now here are some garbage symbols: ?!@#$/\';:[]+=_-.<,>"')

        count = 1
        max_tries = 5
        output = 'first_attempt'
        errors = ['first_attempt', 'slurm_load_jobs', 'slurm_receive_msg', 'send/recv']  # error messages to retry on
        while True in [error in output for error in errors] and count <= max_tries:
            if not output == 'first_attempt':
                time.sleep(30)      # wait 30 seconds before trying again
                count += 1
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read().decode()

        # Use a regular expression to extract the jobid from this string
        pattern = re.compile('[0-9]+')  # todo: it's not inconceivable that this should fail in some cases. Consider moving building this pattern to a method of BatchSystem.
        try:
            return re.findall(pattern, output)[0]
        except IndexError:  # no number in the output
            raise RuntimeError('unable to submit batch job: ' + filename + '\nMessage from batch system: ' + output)

    def cancel_job(self, jobid, settings):
        command = 'scancel ' + str(jobid)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()     # decode converts from bytes-like to string
        return output

    def get_file_suffix(self):
        return 'slurm'


class AdaptPBS(BatchSystem):
    """
    Adapter class for PBS/Torque BatchSystem.

    """

    def get_status(self, jobid, settings):
        if settings.DEBUG:
            return 'C'

        # Iterate through each status three times for robustness in case the status changes while being checked
        for status in ['CE', 'QHW', 'R', 'CE', 'QHW', 'R', 'CE', 'QHW', 'R']:
            command = 'qselect -u $USER -s ' + status
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

            if str(jobid) in output:
                if status == 'CE':      # 'Complete' or 'Exiting'
                    return 'C'
                elif status == 'QHW':   # 'Queued', 'Waiting', or 'Hold'
                    return 'Q'
                elif status == 'R':     # 'Running'
                    return 'R'

        # If this code is reached, it means the jobid did not appear in any of the above searches
        process = subprocess.Popen('qstat -u $USER', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        output = process.stdout.read().decode()  # decode converts from bytes-like to string
        raise RuntimeError('unexpected PBS status for jobid: ' + str(jobid) + '\nFull status string: ' + output)

    def submit_batch(self, filename, settings):
        command = 'qsub {file}'.replace('{file}', filename)

        if settings.DEBUG:
            command = ('echo "this is a nonsense string for testing purposes: 123456, '
                       'now here are some garbage symbols: ?!@#$/\';:[]+=_-.<,>"')

        count = 1
        max_tries = 5
        output = 'first_attempt'
        errors = ['first_attempt', 'send/recv']  # error messages to retry on
        while True in [error in output for error in errors] and count <= max_tries:
            if not output == 'first_attempt':
                time.sleep(30)      # wait 30 seconds before trying again
                count += 1
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
            output = process.stdout.read().decode()

        # Use a regular expression to extract the jobid from this string
        pattern = re.compile('[0-9]+')  # todo: it's not inconceivable that this should fail in some cases. Consider moving building this pattern to a method of BatchSystem.
        try:
            return re.findall(pattern, output)[0]
        except IndexError:  # no number in the output
            raise RuntimeError('unable to submit batch job: ' + filename + '\nMessage from batch system: ' + output)

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

    def get_file_suffix(self):
        return 'pbs'


class AdaptNone(BatchSystem):
    """
    No batch system; jobs run on Unix directly.

    """

    def get_status(self, jobid, settings):
        # 'Q'ueued is impossible; either it's running or it's not
        try:
            proc = psutil.Process(jobid).status()
            if proc in [psutil.STATUS_RUNNING, psutil.STATUS_SLEEPING, psutil.STATUS_DISK_SLEEP]:
                return 'R'
            elif proc in [psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD]:
                return 'C'
            else:
                raise RuntimeError('unexpected process state for job ' + str(jobid) + ': ' + proc)
        except (psutil.NoSuchProcess, ProcessLookupError):
            return 'C'

    def submit_batch(self, filename, settings):
        command = 'bash {file} &'.replace('{file}', filename)

        if settings.DEBUG:
            command = ('echo "this is a nonsense string for testing purposes: 123456, '
                       'now here are some garbage symbols: ?!@#$/\';:[]+=_-.<,>"')

        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       close_fds=True, shell=True)
        pid = process.pid   # use pid in place of batch system jobid

        return pid

    def cancel_job(self, jobid, settings):
        try:
            p = psutil.Process(jobid)   # jobid is actually a pid
            p.kill()
        except psutil.NoSuchProcess:
            pass

    def get_file_suffix(self):
        return 'sh'
