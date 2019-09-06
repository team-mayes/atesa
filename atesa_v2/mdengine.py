"""
Interface for MDEngine objects. New MDEngines can be implemented by constructing a new class that inherits from MDEngine
and implements its abstract methods.
"""

import abc
import pytraj
import os

class MDEngine(abc.ABC):
    """
    Abstract base class for molecular dynamics engines.

    Implements methods for all of the engine-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_last_frame(self, trajectory, settings):
        """
        Return a new file containing just the last frame of a trajectory in Amber .rst7 format

        Parameters
        ----------
        trajectory : str
            Name of trajectory file to obtain last frame from
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        last_frame : str
            Name of .rst7 format coordinate file corresponding to last frame of trajectory

        """
        pass


class AdaptAmber(MDEngine):
    """
    Adapter class for Amber MDEngine.

    """

    def get_last_frame(self, trajectory, settings):
        new_restart_name = trajectory + '_last_frame.rst7'
        traj = pytraj.iterload(trajectory, settings.topology)
        pytraj.write_traj(new_restart_name, traj, format='rst7', frame_indices=[-1], options='multi', overwrite=True)
        try:
            os.rename(new_restart_name + '.1', new_restart_name)
        except OSError:
            if not os.path.exists(new_restart_name):
                raise OSError('expected pytraj to write either ' + new_restart_name + ' or ' + new_restart_name + '.1, '
                              'but found neither.')
        return new_restart_name

