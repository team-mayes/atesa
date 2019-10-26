"""
Interface for MDEngine objects. New MDEngines can be implemented by constructing a new class that inherits from MDEngine
and implements its abstract methods.
"""

import abc
import os
import pytraj
import mdtraj

class MDEngine(abc.ABC):
    """
    Abstract base class for molecular dynamics engines.

    Implements methods for all of the engine-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_frame(self, trajectory, frame, settings):
        """
        Return a new file containing just the frame'th frame of a trajectory in Amber .rst7 format

        Parameters
        ----------
        trajectory : str
            Name of trajectory file to obtain last frame from
        frame : int
            Index of frame to return; -1 gives last frame, 0 is invalid
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        last_frame : str
            Name of .rst7 format coordinate file corresponding to desired frame of trajectory, if it exists; an empty
            string otherwise

        """
        pass


class AdaptAmber(MDEngine):
    """
    Adapter class for Amber MDEngine.

    """

    def get_frame(self, trajectory, frame, settings):
        new_restart_name = trajectory + '_frame_' + str(frame) + '.rst7'
        if not os.path.exists(trajectory):
            return ''   # since it's possible to call this before the trajectory file has been initialized
        if frame >= 1:
            shift_frame = frame - 1     # because write_traj is 0-indexed but get_frame is 1-indexed
        elif frame == -1:
            shift_frame = -1
        else:
            raise IndexError('invalid frame index for get_frame: ' + str(frame) + ' (must be >= 1, or exactly -1)')

        # Use mdtraj to check for non-zero trajectory length (pytraj gives an error below if n_frames = 0)
        try:
            traj = mdtraj.load(trajectory, top=settings.topology)
            if traj.n_frames == 0:
                del traj
                return ''
        except ValueError:      # sometimes this is the result of trying to load a trajectory too early
            return ''

        traj = pytraj.iterload(trajectory, settings.topology)
        try:
            pytraj.write_traj(new_restart_name, traj, format='rst7', frame_indices=[shift_frame], options='multi', overwrite=True)
        except ValueError:  # pytraj raises a ValueError if frame index is out of range
            raise IndexError('frame index ' + str(frame) + ' is out of range for trajectory: ' + trajectory)
        try:
            os.rename(new_restart_name + '.1', new_restart_name)
        except OSError:
            if not os.path.exists(new_restart_name):
                raise OSError('expected pytraj to write either ' + new_restart_name + ' or ' + new_restart_name + '.1, '
                              'but found neither.')

        return new_restart_name

