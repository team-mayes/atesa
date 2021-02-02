"""
Interface for MDEngine objects. New MDEngines can be implemented by constructing a new class that inherits from MDEngine
and implements its abstract methods.
"""

import abc
import os
import pytraj
import mdtraj
import mdtraj.formats

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
            Index of frame to return; 1-indexed, -1 gives last frame, 0 is invalid
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        last_frame : str
            Name of .rst7 format coordinate file corresponding to desired frame of trajectory, if it exists; an empty
            string otherwise

        """
        pass

    @abc.abstractmethod
    def write_find_ts_restraint(self, basin, inp_file):
        """
        Return the input file for a restrained simulation that performs barrier crossing in find_ts

        Parameters
        ----------
        basin : tuple
            Either settings.commit_fwd or settings.commit_bwd, defining the basin to restrain towards
        inp_file : str
            Name of original input file without restraint

        Returns
        -------
        input_file : str
            Name of the simulation input file

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

        # Check for non-zero trajectory length
        try:
            with pytraj.utils.context.capture_stdout():     # suppress C++ errors if this file is empty
                traj = pytraj.iterload(trajectory, settings.topology)
            if traj.n_frames == 0:
                del traj
                return ''
        except ValueError:      # sometimes this is the result of trying to load a trajectory too early
            return ''

        try:
            pytraj.write_traj(new_restart_name, traj, format='rst7', frame_indices=[shift_frame], options='multi', overwrite=True, velocity=True)
        except ValueError:  # pytraj raises a ValueError if frame index is out of range
            raise IndexError('frame index ' + str(frame) + ' is out of range for trajectory: ' + trajectory)
        except AssertionError:  # sometimes there's an assertion error when shift_frame = -1; cause unknown, but this fixes it
            if shift_frame == -1:
                shift_frame = traj.n_frames - 1
            try:
                pytraj.write_traj(new_restart_name, traj, format='rst7', frame_indices=[shift_frame], options='multi',
                                  overwrite=True, velocity=True)
            except ValueError:  # pytraj raises a ValueError if frame index is out of range
                raise IndexError('frame index ' + str(frame) + ' is out of range for trajectory: ' + trajectory)
        try:
            os.rename(new_restart_name + '.1', new_restart_name)
        except OSError:
            if not os.path.exists(new_restart_name):
                raise OSError('expected pytraj to write either ' + new_restart_name + ' or ' + new_restart_name + '.1, '
                              'but found neither.\nThe most likely explanation is that the ATESA process is unable to '
                              'write to the working directory (' + settings.working_directory + '). You may have run '
                              'out of storage space.')

        return new_restart_name

    def write_find_ts_restraint(self, basin, inp_file):
        # Writing an amber-style restraint file from scratch
        open('find_ts_restraints.disang', 'w').write('')  # initialize file
        with open('find_ts_restraints.disang', 'a') as f:
            f.write('DISANG restraint file produced by ATESA with find_ts > 0 to produce transition state guesses\n')
            for def_index in range(len(basin[0])):
                extra = 0  # additional distance to add to basin definition to push *into* basin rather than to its edge
                if basin[3][def_index] == 'lt':
                    extra = -0.1 * basin[2][def_index]  # minus 10% # todo: this doesn't work for negative angles/dihedrals (moves them towards zero instead of "more negative"), which is fine for now since angles and dihedrals are not yet supported
                elif basin[3][def_index] == 'gt':
                    extra = 0.1 * basin[2][def_index]   # plus 10%
                else:
                    raise RuntimeError('entries in the last list in commitment definitions (commit_fwd and commit_bwd) '
                                       'must be either \'lt\' (less than) or \'gt\' (greater than)')

                f.write(' &rst\n')
                f.write('  iat=' + str(basin[0][def_index]) + ',' + str(basin[1][def_index]) + ',\n')
                f.write('  r1=0, r2=' + str(basin[2][def_index] + extra) + \
                        ', r3=' + str(basin[2][def_index] + extra) + \
                        ', r4=' + str(basin[2][def_index] + extra + 2) + ',\n')
                f.write('  rk2=500, rk3=500,\n')
                f.write(' &end\n')

        return inp_file     # unmodified, as find_ts_amber already includes the reference to find_ts_restraints.disang
