"""
Interface for TaskManager objects. New TaskManagers can be implemented by constructing a new class that inherits from
TaskManager and implements its abstract methods.
"""

import abc
try:
    import factory
except ModuleNotFoundError:
    import atesa_v2.factory as factory

class TaskManager(abc.ABC):
    """
    Abstract base class for task managers.

    Implements methods for all of the task manager-specific tasks that ATESA might need.

    """

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


class AdaptSimple(TaskManager):
    """
    Adapter class for my simple, no-frills task manager.

    Just interfaces directly with the batch system through the terminal.

    """

    def submit_batch(self, filename, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        command = batchsystem.get_submit_command().replace('{file}', filename)

        if settings.DEBUG:
            return '123456'

        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   close_fds=True, shell=True)
        return process.stdout.read()
