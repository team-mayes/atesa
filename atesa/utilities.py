"""
Utility functions that don't properly fit anywhere else. These are operations that aren't specific to any particular
interface or script.
"""

import pytraj
import mdtraj
import shutil
import fileinput
import re
import os
import sys
import glob
import math
import numpy
import numpy as np  # to support custom CVs that use this common alias
import pickle
import copy
import subprocess
import warnings
import itertools
from atesa import factory
from multiprocess import Pool
from statsmodels.tsa import stattools

def check_commit(filename, settings, all_frames=False):
    """
    Check commitment of coordinate file to basins defined by settings.commit_fwd and settings.commit_bwd.

    Raises a RuntimeError if the supplied file is not suitable for checking for commitment.

    Parameters
    ----------
    filename : str
        Name of coordinate or trajectory file to be checked. If the file has more than one frame, only the last frame
        will be loaded (unless all_frames is True)
    settings : argparse.Namespace
        Settings namespace object
    all_frames : bool
        If True, returns a list of commit_flags (one for each frame of the trajectory) instead of just one for the last
        frame.

    Returns
    -------
    commit_flag : str
        Either 'fwd' or 'bwd' if the coordinates are in the corresponding basin, or '' if in neither. If all_frames is
        True, then a list of such strings is returned instead.

    """
    # For some really weird reason I need to import these here for them to correctly be available to eval later on
    import mdtraj
    import pytraj
    import numpy

    # Need check_commit to be able to check the status of CP2K trajectories before they've been renamed
    # todo: come up with something better than this kludge
    if settings.md_engine.lower() == 'cp2k':
        if not os.path.exists(filename):
            files = glob.glob('*-' + filename + '-*')
            if len(files) == 1:
                filename = files[0]
            else:
                raise RuntimeError('Could not identify file: ' + filename + '\nIt does not exist and there is no single'
                                   ' unique file matching the pattern \'*-' + filename + '-*\'. This probably indicates'
                                   ' a problem with the simulation that should have created this file.')

    try:
        full_mtraj = mdtraj.load(filename, top=settings.topology)
        if full_mtraj.n_frames == 0:
            raise RuntimeError('Attempted to check commitment for file: ' + filename + ' but it has no frames.')
    except ValueError as e:
        raise RuntimeError('Unable to load file: ' + filename + ' for checking commitment with the specified topology '
                           'file.\nmdtraj returned error: ' + str(e))
    commit_flag = []        # initialize
    need_pytraj = False     # initialize

    if not all_frames:
        full_mtraj = full_mtraj[-1]

    for ii in range(full_mtraj.n_frames):
        mtraj = full_mtraj[ii]
        commit_flag.append('')  # initialize commit_flag for this frame

        # At this point, the code bifurcates: a legacy branch for the nested list interatomic distance definitions
        # formalism and a new branch for the more flexible boolean evaluation formalism. In order to support mixed
        # formalisms, the bifurcation is done twice: once for commit_fwd and once for commit_bwd

        # Legacy fwd
        if len(settings.commit_fwd) == 4 and type(settings.commit_fwd[0]) == list:
            for i in range(len(settings.commit_fwd[2])):
                if settings.commit_fwd[3][i] == 'lt':
                    if mdtraj.compute_distances(mtraj, numpy.array([[settings.commit_fwd[0][i] - 1, settings.commit_fwd[1][i] - 1]]))[0][0] * 10 <= settings.commit_fwd[2][i]:
                        commit_flag[-1] = 'fwd'     # if a committor test is passed, testing moves on to the next one.
                    else:
                        commit_flag[-1] = ''
                        break                   # if a committor test is not passed, all testing in this direction fails
                elif settings.commit_fwd[3][i] == 'gt':
                    if mdtraj.compute_distances(mtraj, numpy.array([[settings.commit_fwd[0][i] - 1, settings.commit_fwd[1][i] - 1]]))[0][0] * 10 >= settings.commit_fwd[2][i]:
                        commit_flag[-1] = 'fwd'
                    else:
                        commit_flag[-1] = ''
                        break
                else:
                    raise ValueError('An incorrect committor definition \"' + settings.commit_fwd[3][i] + '\" was given for '
                                     'index ' + str(i) + ' in the \'fwd\' direction.')

        # New fwd
        else:
            # Need to set variables for evaluation: pytraj traj, and traj_name (mdtraj mtraj is already loaded)
            # First, evaluate boolean determining whether we need to bother loading pytraj
            need_pytraj = any([('traj' in cv and not 'mtraj' in cv) for cv in settings.commit_fwd])
            if need_pytraj:
                traj = pytraj.load(filename, settings.topology)
            traj_name = filename
            scope = locals()
            if all([eval(item, scope) for item in settings.commit_fwd]):
                commit_flag[-1] = 'fwd'
            else:
                commit_flag[-1] = ''

        if commit_flag[-1] == '':  # only bother checking for bwd commitment if not fwd committed
            # Legacy bwd
            if len(settings.commit_bwd) == 4 and type(settings.commit_bwd[0]) == list:
                for i in range(len(settings.commit_bwd[2])):
                    if settings.commit_bwd[3][i] == 'lt':
                        if mdtraj.compute_distances(mtraj, numpy.array([[settings.commit_bwd[0][i] - 1, settings.commit_bwd[1][i] - 1]]))[0][0] * 10 <= settings.commit_bwd[2][i]:
                            commit_flag[-1] = 'bwd' # if a committor test is passed, testing moves on to the next one.
                        else:
                            commit_flag[-1] = ''
                            break               # if a committor test is not passed, all testing in this direction fails
                    elif settings.commit_bwd[3][i] == 'gt':
                        if mdtraj.compute_distances(mtraj, numpy.array([[settings.commit_bwd[0][i] - 1, settings.commit_bwd[1][i] - 1]]))[0][0] * 10 >= settings.commit_bwd[2][i]:
                            commit_flag[-1] = 'bwd'
                        else:
                            commit_flag[-1] = ''
                            break
                    else:
                        raise ValueError('An incorrect committor definition \"' + settings.commit_bwd[3][i] + '\" was given for'
                                         ' index ' + str(i) + ' in the \'bwd\' direction.')

            # New bwd
            else:
                # Need to set variables for evaluation: pytraj traj, and traj_name (mdtraj mtraj is already loaded)
                # First, evaluate boolean determining whether we need to bother loading pytraj
                if not need_pytraj:     # skip this if need_pytraj has already evaluated as true
                    need_pytraj = any([('traj' in cv and not 'mtraj' in cv) for cv in settings.commit_bwd])
                    if need_pytraj:
                        traj = pytraj.load(filename, settings.topology)
                traj_name = filename
                scope = locals()
                if all([eval(item, scope) for item in settings.commit_bwd]):
                    commit_flag[-1] = 'bwd'
                else:
                    commit_flag[-1] = ''

    if not all_frames:
        assert len(commit_flag) == 1
        commit_flag = commit_flag[0]

    return commit_flag


def get_cvs(filename, settings, reduce=False, frame=0, only_in_rc=False):
    """
    Get CV values for a coordinate file given by filename, as well as rates of change if settings.include_qdot = True.

    If reduce = True, the returned CVs will be reduced to between 0 and 1 based on the minimum and maximum values of
    that CV in as.out, which is assumed to exist at settings.as_out_file.

    If frame is > 0 or 'all', filename will be interpreted as a trajectory file instead of a .rst7 coordinate file.
    This option is incompatible with settings.include_qdot = True.

    Parameters
    ----------
    filename : str
        Name of .rst7-formatted coordinate file to be checked (or trajectory file if frame > 0 or frame == 'all')
    settings : argparse.Namespace
        Settings namespace object
    reduce : bool
        Boolean determining whether to reduce the CV values for use in evaluating an RC value that uses reduced values
    frame : int or str
        Frame of trajectory filename to use (0-indexed, negatives read from end); or, 'all' returns the CVs for each
        frame of the input trajectory, separated by newlines. 'all' is incompatible with settings.include_qdot = True
    only_in_rc : bool
        If this is True, then do not evaluate and return 0 for all CVs not present in the RC defined by
        settings.rc_definition. This is useful to speed up evaluation when all CVs are not necessary. Applies also to
        corresponding qdot values as appropriate.

    Returns
    -------
    output : str
        Space-separated list of CV values for the given coordinate file

    """

    def increment_coords():
        # Produces a new coordinate file from the existing one by incrementing coordinates by their rates of change
        # Returns the name of the newly-created coordinate file
        byline = open(filename).readlines()
        pattern = re.compile('-*[0-9.]+')           # regex to match numbers including decimals and negatives
        n_atoms = pattern.findall(byline[1])[0]     # number of atoms indicated on second line of .rst file

        shutil.copyfile(filename, filename + '_temp.rst7')
        for i, line in enumerate(fileinput.input(filename + '_temp.rst7', inplace=1)):
            if int(n_atoms)/2 + 2 > i >= 2:
                newline = line
                coords = pattern.findall(newline)                                          # line of coordinates
                try:
                    vels = pattern.findall(byline[i + int(math.ceil(int(n_atoms)/2))])     # corresponding velocities
                except IndexError:
                    os.remove(filename + '_temp.rst7')      # to clean up
                    os.remove(filename + '_temp.rst7.bak')  # to clean up
                    fileinput.close()
                    raise IndexError('get_cvs.increment_coords() encountered an IndexError. This is caused '
                             'by attempting to read qdot values from a coordinate file lacking velocity information, or'
                             ' else by that file being truncated. Ensure that the relevant simulation input file is set'
                             ' to write velocities to output trajectories, as this may not be default behavior. In '
                             'Amber, this is accomplished by setting ntwv=-1 in the MD input file. Alternatively, if '
                             'you don\'t want to include velocity terms in the output file, you can set include_qdot = '
                             'False in the ATESA config file.\n'
                             'The offending file is: ' + filename)

                # Sometimes items in coords or vels 'stick together' at a negative sign (e.g., '-1.8091748-112.6420521')
                # This next loop is just to split them up
                for index in range(len(coords)):
                    length = len(coords[index])                     # length of string representing this coordinate
                    replace_string = str(float(coords[index]) + float(vels[index]))[0:length-1]
                    while len(replace_string) < length:
                        replace_string += '0'
                    newline = newline.replace(coords[index], replace_string)
                sys.stdout.write(newline)
            else:
                sys.stdout.write(line)

        return filename + '_temp.rst7'

    def reduce_cv(unreduced_value, local_index, rc_minmax):
        # Returns a reduced value for a CV given an unreduced value and the index within as.out corresponding to that CV
        this_min = rc_minmax[0][local_index]
        this_max = rc_minmax[1][local_index]
        if this_min == this_max:
            return unreduced_value      # can't reduce when a CV only has a single value across rc_minmax
        return (float(unreduced_value) - this_min) / (this_max - this_min)

    def cv_is_in_rc(cv_index):
        # Helper function checks if the zero-indexed cv_index corresponds to a CV in settings.rc_definition
        if settings.rc_definition == '':
            raise RuntimeError('utilities.get_cvs was called with only_in_rc=True, but rc_definition was empty. You '
                               'must provide the RC definition in the configuration file to use this option.')
        # Have to be careful here not to interpret, e.g., "CV10" as matching "CV1", so we need to be a bit clever
        pattern = re.compile('CV[0-9]+')
        cvs_in_rc = pattern.findall(settings.rc_definition)
        cv_to_look_for = 'CV' + str(int(cv_index) + 1)   # cv_index is zero-indexed, but CVs in RCs are 1-indexed
        return cv_to_look_for in cvs_in_rc

    if frame == 'all' and settings.include_qdot:
        raise RuntimeError('utilities.get_cvs cannot be called with frame=all if settings.include_qdot == True')

    if not os.getcwd() == settings.working_directory:
        os.chdir(settings.working_directory)    # make sure we're in the working directory

    # boolean determining whether we need to bother loading pytraj
    need_pytraj = any([('traj' in cv and not 'mtraj' in cv) for cv in settings.cvs])

    # define variables as they are expected to be interpreted in CV definitions (these are potentially invoked by eval)
    delete_traj_name = False  # so we don't delete the traj_name file later unless we created it below
    if frame in [0, 'all']:
        try:
            full_mtraj = mdtraj.load(filename, top=settings.topology)
        except TypeError:   # can happen during resampling if trajectory finished after controlling ATESA job did
            mdengine = factory.mdengine_factory(settings.md_engine)
            mdengine.control_rst_format(filename)
            full_mtraj = mdtraj.load(filename, top=settings.topology)
        if need_pytraj:  # only load pytraj if it's needed
            full_traj = pytraj.iterload(filename, settings.topology)
        traj_name = filename
    else:
        if need_pytraj:
            if frame < 0:
                full_traj = pytraj.iterload(filename, settings.topology)
                pframe = full_traj.n_frames + frame
                full_traj = pytraj.iterload(filename, settings.topology, frame_slice=[(pframe, pframe + 1)])
            else:
                full_traj = pytraj.iterload(filename, settings.topology, frame_slice=[(frame, frame + 1)])
        full_mtraj = mdtraj.load(filename, top=settings.topology)[frame]
        if True in ['traj_name' in cv for cv in settings.cvs]:  # if True, need .rst7 formatted file to operate on
            mdengine = factory.mdengine_factory(settings.md_engine)
            traj_name = mdengine.get_frame(filename, frame, settings)
            delete_traj_name = True     # to be sure we clean up traj_name later

    rc_minmax = [[],[]]
    if reduce:
        # Try to load from file if available
        size = os.path.getsize(settings.as_out_file)
        if os.path.exists(str(size) + '_minmax.pkl'):
            rc_minmax = pickle.load(open(str(size) + '_minmax.pkl', 'rb'))
        else:
            # Prepare new cv_minmax list
            asout_lines = [[float(item) for item in line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')] for line in open(settings.as_out_file, 'r').readlines()]
            open(settings.as_out_file, 'r').close()
            mapped = list(map(list, zip(*asout_lines)))
            rc_minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]
            pickle.dump(rc_minmax, open(str(size) + '_minmax.pkl', 'wb'), protocol=2)

    output = ''
    values = []
    sigfigs = '%.' + str(settings.sigfigs) + 'E'

    if frame == 'all':
        n_iter = full_mtraj.n_frames
    else:
        n_iter = 1

    for iter in range(n_iter):  # iterate over all frames
        local_index = -1
        # set traj and mtraj to only the desired frame; has no (important) effect if there's only one frame anyway
        mtraj = full_mtraj[iter]
        if need_pytraj:     # only load traj if it's needed
            traj = full_traj[iter:iter + 1:1]
        for cv in settings.cvs:
            local_index += 1
            if (not only_in_rc) or cv_is_in_rc(local_index):
                evaluation = eval(cv)       # evaluate the CV definition code (potentially using traj, mtraj, and/or traj_name)
                if settings.include_qdot:   # want to save values for later
                    values.append(float(evaluation))
                if reduce:                  # reduce values if necessary
                    evaluation = reduce_cv(evaluation, local_index, rc_minmax)
            else:
                evaluation = 0
            output += sigfigs % float(evaluation) + ' '
        if n_iter > 1:
            output += '\n'
        elif settings.include_qdot and not settings.job_type == 'equilibrium_path_sampling':  # if True, then we want to include rate of change for every CV, too
            # Strategy here is to write a new temporary .rst7 file by incrementing all the coordinate values by their
            # corresponding velocity values, load it as a new iterload object, and then rerun our analysis on that.
            incremented_filename = increment_coords()
            temp_traj_name = traj_name
            traj_name = incremented_filename    # to make available to cv evaluation code as needed
            mtraj = mdtraj.load(incremented_filename, top=settings.topology)
            if need_pytraj:  # only load traj if it's needed
                traj = pytraj.iterload(incremented_filename, settings.topology)
            local_index = -1
            for cv in settings.cvs:
                local_index += 1
                if (not only_in_rc) or cv_is_in_rc(local_index):
                    evaluation = eval(cv) - values[local_index]  # Subtract value 1/20.455 ps earlier from value of cv
                    if reduce:
                        try:
                            evaluation = reduce_cv(evaluation, local_index + len(settings.cvs), rc_minmax)
                        except IndexError:
                            raise RuntimeError('attempted to obtain a reduced value for the rate-of-change ("qdot") of CV' +
                                               str(local_index) + ', but the corresponding column appears to be missing from '
                                               'the aimless shooting output file: ' + settings.as_out_file + '. If this is the '
                                               'wrong file, please specify the right one with the "as_out_file" option in the '
                                               'configuration file. If this is the right file but it does not contain qdot '
                                               'terms, please add "include_qdot = False" to the configuration file.')
                else:
                    evaluation = 0
                output += sigfigs % float(evaluation) + ' '
            os.remove(incremented_filename)     # clean up temporary file
            traj_name = temp_traj_name

    while output[-1] in [' ', '\n']:
        output = output[:-1]                # remove trailing space and/or newline on terminating line
    output = output.replace(' \n', '\n')    # remove trailing spaces on non-terminating lines

    if delete_traj_name:
        os.remove(traj_name)

    return output


def rev_vels(restart_file):
    """
    Reverse all the velocity terms in a restart file and return the name of the new, 'reversed' file.

    Parameters
    ----------
    restart_file : str
        Filename of the 'fwd' restart file, in .rst7 format

    Returns
    -------
    reversed_file : str
        Filename of the newly written 'bwd' restart file, in .rst7 format

    """

    byline = open(restart_file).readlines()
    open(restart_file).close()
    pattern = re.compile(r'[-0-9.]+')            # regex to match numbers including decimals and negatives
    pattern2 = re.compile(r'\s[-0-9.]+')         # regex to match numbers including decimals and negatives, with one space in front
    try:
        n_atoms = pattern.findall(byline[1])[0]     # number of atoms indicated on second line of .rst7 file
    except IndexError:
        raise RuntimeError('restart file: ' + restart_file + ' exists but is improperly formatted. It must contain the '
                           'number of atoms as the first integer on its second line.')
    offset = 2                  # appropriate for n_atoms is odd; offset helps avoid modifying the box line
    if int(n_atoms) % 2 == 0:   # if n_atoms is even...
        offset = 1              # appropriate for n_atoms is even

    try:
        name = restart_file[:restart_file.rindex('.')]  # everything before last '.', to remove file extension
    except ValueError:
        name = restart_file     # if no '.' in the filename

    shutil.copyfile(restart_file, name + '_bwd.rst7')
    for i, line in enumerate(fileinput.input(name + '_bwd.rst7', inplace=1)):
        if int(n_atoms) / 2 + 2 <= i <= int(n_atoms) + offset:  # if this line is a velocity line
            newline = line
            for vel in pattern2.findall(newline):
                if '-' in vel:
                    newline = newline.replace(vel, '  ' + vel[2:], 1)   # replace ' -magnitude' with '  magnitude'
                else:
                    newline = newline.replace(vel, '-' + vel[1:], 1)    # replace ' magnitude' with '-magnitude'
            sys.stdout.write(newline)
        else:  # if not a velocity line
            sys.stdout.write(line)

    return name + '_bwd.rst7'


def evaluate_rc(rc_definition, cv_list):
    """
    Evaluate the RC value given by RC definition for the given list of CV values given by cv_list.

    Parameters
    ----------
    rc_definition : str
        A reaction coordinate definition formatted as a string of python-readable code with "CV[X]" standing in for the
        Xth CV value (zero-indexed); e.g., "CV2" has value "4" in the cv_list [1, -2, 4, 6]
    cv_list : list
        A list of CV values whose indices correspond to the desired values in rc_definition

    Returns
    -------
    rc_value : float
        The value of the reaction coordinate given the values in cv_list

    """

    # Fill in CV[X] slots with corresponding values from cv_list
    for i in reversed(range(len(cv_list))):     # reversed so that e.g. CV10 doesn't get interpreted as '[CV1]0'
        rc_definition = rc_definition.replace('CV' + str(i + 1), str(cv_list[i]))

    # Evaluate the filled-in rc_definition and return the result
    return eval(rc_definition)


def resample(settings, partial=False, full_cvs=False, only_full_cvs=False):
    """
    Resample each shooting point in each thread with different CV definitions to produce new output files with extant
    aimless shooting data.

    This function also assesses decorrelation times and produces one or more decorrelated output files. If and only if
    settings.information_error_checking == True, decorrelated files are produced at each settings.information_error_freq
    increment. In this case, if partial == True, decorrelation will only be assessed for data lengths absent from the
    info_err.out file in the working directory.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object
    partial : bool
        If True, reads the info_err.out file and only builds new decorrelated output files where the corresponding lines
        are missing from that file. If partial == False, decorrelation is assessed for every valid data length. Has no
        effect if not settings.information_error_checking.
    full_cvs : bool
        If True, also resamples as_full_cvs.out using every prod trajectory in the working directory.

    Returns
    -------
    None

    """

    # todo: test this more thoroughly using a dummy thread and a manual decorrelation time calculation using different software

    # This function is sometimes called from outside the working directory, so make sure we're there
    os.chdir(settings.working_directory)

    # Load in allthreads from restart.pkl
    try:
        allthreads = pickle.load(open('restart.pkl', 'rb'))
    except FileNotFoundError:
        raise FileNotFoundError('resample = True requires restart.pkl, but could not find one in working directory: '
                                + settings.working_directory)

    if not only_full_cvs:
        # Remove pre-existing output files if any, initialize new one
        open(settings.working_directory + '/as_raw_resample.out', 'w').close()
        if settings.information_error_checking:
            open(settings.working_directory + '/as_raw_timestamped.out', 'w').close()

        # Open files for writing outside loop (much faster than opening/closing for each write)
        f1 = open(settings.working_directory + '/as_raw_resample.out', 'a')
        if settings.information_error_checking:
            f2 = open(settings.working_directory + '/as_raw_timestamped.out', 'a')

        # Iterate through each thread's history.init_coords list and obtain CV values as needed
        for thread in allthreads:
            thread.this_cvs_list = []       # initialize full nested list of CV values for this thread
            thread.cvs_for_later = []       # need this one with empty lists for failed moves, for indexing reasons
            for step_index in range(len(thread.history.prod_results)):
                if thread.history.prod_results[step_index][0] in ['fwd', 'bwd']:
                    if thread.history.prod_results[step_index][0] == 'fwd':
                        this_basin = 'B'
                    else:  # 'bwd'
                        this_basin = 'A'

                    # Get CVs for this shooting point   # todo: a bit sloppy... can I clean this up?
                    try:
                        if not os.path.exists(thread.history.init_coords[step_index][0]):
                            warnings.warn('attempted to resample ' + thread.history.init_coords[step_index][0] + ' but no such '
                                          'file exists in the working directory\nSkipping and continuing', RuntimeWarning)
                            thread.cvs_for_later.append([])
                            continue        # skip to next step_index
                    except IndexError:  # getting cv's failed (maybe corrupt coordinate file) so consider this step failed
                        thread.cvs_for_later.append([])
                        continue        # skip to next step_index
                    try:
                        this_cvs = get_cvs(thread.history.init_coords[step_index][0], settings)
                    except IndexError:  # getting cv's failed (maybe corrupt coordinate file) so consider this step failed
                        thread.cvs_for_later.append([])
                        continue        # skip to next step_index

                    # Write CVs to as_raw_resample.out
                    f1.write(this_basin + ' <- ' + this_cvs + '\n')
                    f1.flush()
                    if settings.information_error_checking:
                        f2.write(str(thread.history.timestamps[step_index]) + ' ' + this_basin + ' <- ' + this_cvs + '\n')
                        f2.flush()

                    # Append this_cvs to running list for evaluating decorrelation time
                    thread.this_cvs_list.append([[float(item) for item in this_cvs.split(' ')], thread.history.timestamps[step_index]])
                    thread.cvs_for_later.append([float(item) for item in this_cvs.split(' ')])
                else:
                    thread.cvs_for_later.append([])

        # Close files just to be sure
        f1.close()
        if settings.information_error_checking:
            f2.close()

        if settings.information_error_checking:   # sort timestamped output file
            shutil.copy(settings.working_directory + '/as_raw_timestamped.out', settings.working_directory + '/as_raw_timestamped_copy.out')
            open(settings.working_directory + '/as_raw_timestamped.out', 'w').close()
            with open(settings.working_directory + '/as_raw_timestamped_copy.out', 'r') as f:
                for line in sorted(f):
                    open(settings.working_directory + '/as_raw_timestamped.out', 'a').write(line)
                open(settings.working_directory + '/as_raw_timestamped.out', 'a').close()
            os.remove(settings.working_directory + '/as_raw_timestamped_copy.out')

        # Construct list of data lengths to perform decorrelation for
        if settings.information_error_checking:
            if not partial:
                lengths = [leng for leng in range(settings.information_error_freq, len(open(settings.working_directory + '/as_raw_timestamped.out', 'r').readlines()) + 1, settings.information_error_freq)]
            else:   # if partial
                lengths = [leng for leng in range(settings.information_error_freq, len(open(settings.working_directory + '/as_raw_timestamped.out', 'r').readlines()) + 1, settings.information_error_freq) if not leng in [int(line.split(' ')[0]) for line in open(settings.working_directory + '/info_err.out', 'r').readlines()]]
            pattern = re.compile('[0-9]+')  # pattern for reading out timestamp from string
        else:
            lengths = [len(open(settings.working_directory + '/as_raw_resample.out', 'r').readlines())]
            pattern = None

        # Assess decorrelation and write as_decorr.out
        for length in lengths:
            if settings.information_error_checking:
                suffix = '_' + str(length)     # only use-case with multiple lengths, so this keeps them from stepping on one another's toes
                cutoff_timestamp = int(pattern.findall(open(settings.working_directory + '/as_raw_timestamped.out', 'r').readlines()[length - 1])[0])
            else:
                cutoff_timestamp = math.inf
                suffix = ''
            open(settings.working_directory + '/as_decorr' + suffix + '.out', 'w').close()
            f3 = open(settings.working_directory + '/as_decorr' + suffix + '.out', 'a')
            for thread in allthreads:
                if thread.this_cvs_list:       # if there were any 'fwd' or 'bwd' results in this thread
                    mapped = list(map(list, zip(*[item[0] for item in thread.this_cvs_list if item[1] <= cutoff_timestamp])))   # list of lists of values of each CV

                    slowest_lag = -1    # initialize running tally of slowest autocorrelation time among CVs in this thread
                    if settings.include_qdot:
                        ndims = len(thread.this_cvs_list[0]) / 2   # number of non-rate-of-change CVs
                        if not ndims % 1 == 0:
                            raise ValueError('include_qdot = True, but an odd number of dimensions were found in the '
                                             'threads in restart.pkl, so they can\'t contain inertial terms.')
                        ndims = int(ndims)
                    else:
                        ndims = len(thread.this_cvs_list[0])

                    for dim_index in range(ndims):
                        slowest_lag = -1
                        if mapped:
                            this_cv = mapped[dim_index]
                        if len(this_cv) > 1:
                            this_autocorr = stattools.acf(this_cv, nlags=len(this_cv) - 1, fft=True)
                            for lag in range(len(this_cv) - 1):
                                corr = this_autocorr[lag]
                                if abs(corr) <= 1.96 / numpy.sqrt(len(this_cv)):
                                    slowest_lag = lag + 1
                                    break

                    if slowest_lag > 0:     # only proceed to writing decorrelated output file if a slowest_lag was found
                        # Write the same way as to as_raw_resample.out above, but starting the range at slowest_lag
                        for step_index in range(slowest_lag, len(thread.history.prod_results)):
                            if thread.history.prod_results[step_index][0] in ['fwd', 'bwd'] and thread.history.timestamps[step_index] <= cutoff_timestamp:
                                if thread.history.prod_results[step_index][0] == 'fwd':
                                    this_basin = 'B'
                                else:  # 'bwd'
                                    this_basin = 'A'

                                # Get CVs for this shooting point and write them to the decorrelated output file
                                if thread.cvs_for_later[step_index]:
                                    this_cvs = thread.cvs_for_later[step_index]    # retrieve CVs from last evaluation
                                    f3.write(this_basin + ' <- ' + ' '.join([str(item) for item in this_cvs]) + '\n')

            f3.close()

        # Move resample raw output file to take its place as the only raw output file
        shutil.move(settings.working_directory + '/as_raw_resample.out', settings.working_directory + '/as_raw.out')
    else:
        if not full_cvs:
            raise RuntimeError('resample called with only_full_cvs = True but full_cvs = False, so this does nothing. '
                               'Change one of them!')

    # Implement full_cvs
    if full_cvs:
        open(settings.working_directory + '/as_full_cvs.out', 'w').close()
        temp_settings = copy.deepcopy(settings)
        temp_settings.include_qdot = False  # never want to include_qdot in this upcoming call to get_cvs
        try:
            affinity = len(os.sched_getaffinity(0))
        except AttributeError:  # os.sched_getaffinity raises AttributeError on non-UNIX systems.
            affinity = 1
        if affinity == 1:
            with open(settings.working_directory + '/as_full_cvs.out', 'a') as f:
                for thread in allthreads:
                    for step_index in range(min([len(thread.history.prod_results), len(thread.history.prod_trajs)])):    # just in case one got an extra write in over the other
                        if thread.history.prod_results[step_index] in [['fwd', 'bwd'], ['bwd', 'fwd']]:     # if step accepted
                            for job_index in range(2):
                                if os.path.exists(thread.history.prod_trajs[step_index][job_index]):
                                    f.write(get_cvs(thread.history.prod_trajs[step_index][job_index], temp_settings, False, 'all') + '\n')
        else:   # affinity > 1
            # Map partial_full_cvs calls to available processes
            with Pool(affinity) as p:
                p.starmap(partial_full_cvs, zip(allthreads, ['partial_full_cvs_' + str(thread_index) + '.out' for thread_index in range(len(allthreads))], itertools.repeat(temp_settings)))
            # Finally, combine the partial files into the full file
            with open(settings.working_directory + '/as_full_cvs.out', 'w') as outfile:
                for fname in ['partial_full_cvs_' + str(thread_index) + '.out' for thread_index in range(len(allthreads))]:
                    with open(fname) as infile:
                        for line in infile:
                            if line:    # skip blank lines
                                outfile.write(line)
                    os.remove(fname)

def partial_full_cvs(thread, filename, settings):
    """
    Write the full CVs list to the file given by filename for each accepted move in each thread in threads.

    A helper function for parallelizing get_cvs (has to be defined at top-level for multiprocessing)

    Parameters
    ----------
    thread : Thread
        Thread object to evaluate
    filename : str
        Name of the file to write CVs to
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """
    open(settings.working_directory + '/' + filename, 'w').close()
    with open(settings.working_directory + '/' + filename, 'a') as f:
        for step_index in range(min([len(thread.history.prod_results), len(thread.history.prod_trajs)])):  # just in case one got an extra write in over the other
            if thread.history.prod_results[step_index] in [['fwd', 'bwd'], ['bwd', 'fwd']]:  # if step accepted
                for job_index in range(2):
                    if os.path.exists(thread.history.prod_trajs[step_index][job_index]):
                        f.write(get_cvs(thread.history.prod_trajs[step_index][job_index], settings, False, 'all') + '\n')

def interpret_cv(cv_index, settings):
    """
    Read the cv_index'th CV from settings.cvs, identify its type (distance, angle, dihedral, or differance-of-distances)
    and the atom indices the define it (one-indexed) and return these.

    This function is designed for use in the umbrella_sampling jobtype only. For this reason, it only supports the
    aforementioned CV types. If none of these types appears to fit, this function raises a RuntimeError.

    Parameters
    ----------
    cv_index : int
        The index for the CV to use; e.g., 6 corresponds to CV6. Must be in the range [1, len(settings.cvs)].
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    atoms : list
        A list of 1-indexed atom indices as strings that define the given CV
    optype : str
        A string (either 'distance', 'angle', 'dihedral', or 'diffdistance') corresponding to the type for this CV.
    nat : int
        The number of atoms constituting the given CV

    """
    if not 1 <= cv_index <= len(settings.cvs):
        raise RuntimeError('called interpret_cv with an index outside the range [1, len(settings.cvs)].\n'
                           'len(settings.cvs) = ' + str(len(settings.cvs)) + ' and cv_index = ' + str(cv_index))

    this_cv = settings.cvs[cv_index - 1]  # -1 because CVs are 1-indexed
    if 'pytraj.dihedral' in this_cv or 'mdtraj.compute_dihedrals' in this_cv:
        optype = 'dihedral'
        nat = 4
    elif 'pytraj.angle' in this_cv or 'mdtraj.compute_angles' in this_cv:
        optype = 'angle'
        nat = 3
    elif 'pytraj.distance' in this_cv or 'mdtraj.compute_distances' in this_cv:
        if '-' in this_cv and (this_cv.count('pytraj.distance') == 2 or this_cv.count('mdtraj.compute_distances') == 2 or
                                   ('mdtraj.distance' in this_cv and 'pytraj.distance' in this_cv)):
            optype = 'diffdistance'
            nat = 4
        else:
            optype = 'distance'
            nat = 2
    else:
        raise RuntimeError('unable to discern CV type for CV' + str(int(cv_index)) + '\nOnly '
                           'distances, angles, dihedrals, and differences of distances (all defined '
                           'using either pytraj or mdtraj distance, angle, and/or dihedral functions) '
                           'are supported in umbrella sampling. The offending CV is defined as: ' +
                           this_cv)

    # Get atom indices as a string, then convert to list
    atoms = ''
    if not optype == 'diffdistance':
        count = 0
        for match in re.finditer('[\[\\\']([@0-9]+[,\ ]){' + str(nat - 1) + '}[@0-9]+[\]\\\']',
                                 this_cv.replace(', ', ',')):
            atoms += this_cv.replace(', ', ',')[match.start():match.end()]  # should be only one match
            count += 1
        if not count == 1:
            raise RuntimeError('failed to identify atoms constituting CV definition: ' + this_cv +
                               '\nInterpreted as a ' + optype + ' but found ' + str(count) +
                               ' blocks of atom indices with length ' + str(nat) + '(should be one). Is'
                                                                                   ' this CV formatted in an unusual way?')
        if not atoms:
            raise RuntimeError('unable to identify atoms constituting CV definition: ' + this_cv +
                               '\nIs it formatted in an unusual way?')
        atoms = atoms.replace('[', '').replace(']', ',').replace('\'',
                                                                 '')  # included delimeters for safety, but don't want them
        if '@' in atoms:
            atoms = [item.replace('@', '') for item in atoms.split(' @')]  # pytraj style atom indices
        else:
            atoms = atoms.split(',')  # mdtraj style atom indices
            atoms = [str(int(item) + 1) for item in atoms if not item == '']  # fix zero-indexing in mdtraj
    else:
        count = 0
        for match in re.finditer('[\[\\\']([@0-9]+[,\ ]){1}[@0-9]+[\]\\\']', this_cv.replace(', ', ',')):
            atoms += this_cv.replace(', ', ',')[match.start():match.end()]  # should be two matches
            count += 1
        if not count == 2:
            raise RuntimeError('failed to identify atoms constituting CV definition: ' + this_cv +
                               '\nInterpreted as a difference of distances but found ' + str(count) +
                               ' blocks of atom indices with length 2 (should be two). Is this CV '
                               'formatted in an unusual way?')
        if not atoms:
            raise RuntimeError('unable to identify atoms constituting CV definition: ' + this_cv +
                               '\nIs it formatted in an unusual way?')
        atoms = atoms.replace('[', '').replace(']', ',').replace('\'', '')  # included delimeters for safety, but don't want them
        if '@' in atoms:
            atoms = [item.replace('@', '') for item in atoms.replace('\'', ' ').replace('@', ' @').split()]
        else:
            atoms = atoms.split(',')
            atoms = [str(int(item) + 1) for item in atoms if not item == '']  # fix zero-indexing in mdtraj

    while '' in atoms:
        atoms.remove('')  # remove empty list elements if present

    return atoms, optype, nat


def plot_cvs_vs_rc(full_cvs_file, rc_definition, as_output_file):
    """
    Produce a plot of the values of each CV in the RC as a function of RC value from either aimless shooting or umbrella
    sampling data.

    The plot is made in matplotlib but also printed as text to an output file.

    Parameters
    ----------
    full_cvs_file : str
        Path to full_cvs output file (produced by utilities.resample or resampling.resample_umbrella_sampling)
    rc_definition : str
        String representing the definition of the reduced reaction coordinate
    as_output_file : str
        Path to the as_output_file file passed to lmax.py when generating the reaction coordinate

    Returns
    -------
    None

    """
    # Identify relevant CVs
    pattern = re.compile('CV[0-9]+')
    cv_names = pattern.findall(rc_definition)    # get list of CV names, e.g., ['CV6', 'CV35', ...]
    cvs = [int(cv[2:]) - 1 for cv in cv_names]   # get zero-indexed CV indices

    # Get values of relevant CVs from each line in full_cvs_file
    values = [[] for _ in range(len(cvs))]
    for line in open(full_cvs_file, 'r').readlines():
        ii = -1
        for cv in cvs:
            ii += 1
            values[ii].append(float(line.split()[cv]))     # [[obs1cv1, obs2cv1], [obs1cv2, obs2cv2]]

    # Get data from input file, and determine minimum and maximum values for each CV, get reduced CV values
    input_data = [[float(item) for item in
                   line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')]
                  for line in open(as_output_file, 'r').readlines()]       # [[obs1cv1, obs1cv2], [obs2cv1, obs2cv2]]
    mapped = list(map(list, zip(*input_data)))                  # [[obs1cv1, obs2cv1], [obs1cv2, obs2cv2]]
    minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]    # [[mincv1, mincv2], [maxcv1, maxcv2]]
    n_cvs = len(minmax[0])  # number of CVs total in as_output_file, for use later
    minmax = [[minmax[0][ii] for ii in cvs], [minmax[1][ii] for ii in cvs]]     # cut down to just the relevant values
    reduced_values = [[(values[ii][jj] - minmax[0][ii]) / (minmax[1][ii] - minmax[0][ii]) for ii in range(len(values))] for jj in range(len(values[0]))]
    reduced_values = list(map(list, zip(*reduced_values)))  # transpose to match shape of values

    rcs = []     # list of RC values for each sample
    for jj in range(len(reduced_values[0])):    # iterate over samples
        cvs_list = [0 for ii in range(n_cvs)]    # evaluate_rc expects an entry for every CV, not just those in the RC
        for kk in range(len(cvs)):
            cvs_list[cvs[kk]] = reduced_values[kk][jj]
        rcs.append(evaluate_rc(rc_definition, cvs_list))

    # Bin each RC value and evaluate average and stdev of each CV value in each bin
    bins = numpy.linspace(min(rcs), max(rcs), 20) # todo: should the number of bins be a setting?

    assert len(values[0]) == len(rcs)   # sanity check

    means = [[numpy.mean([values[cv][jj] for jj in range(len(rcs)) if bins[bin_ii] <= rcs[jj] <= bins[bin_ii + 1]]) for
              bin_ii in range(len(bins) - 1)] for cv in range(len(values))]
    stds = [[numpy.std([values[cv][jj] for jj in range(len(rcs)) if bins[bin_ii] <= rcs[jj] <= bins[bin_ii + 1]]) for
              bin_ii in range(len(bins) - 1)] for cv in range(len(values))]
    reduced_means = [[numpy.mean([reduced_values[cv][jj] for jj in range(len(rcs)) if bins[bin_ii] <= rcs[jj] <= bins[bin_ii + 1]]) for
              bin_ii in range(len(bins) - 1)] for cv in range(len(values))]
    reduced_stds = [[numpy.std([reduced_values[cv][jj] for jj in range(len(rcs)) if bins[bin_ii] <= rcs[jj] <= bins[bin_ii + 1]]) for
              bin_ii in range(len(bins) - 1)] for cv in range(len(values))]

    # Write results to output file and plot in both reduced and non-reduced forms
    with open('cvs_vs_rc.out', 'w') as f:
        f.write('RC value')
        for name in cv_names:
            f.write('; ' + name + ' mean')
            f.write('; ' + name + ' st.dev.')
        for name in cv_names:
            f.write('; reduced ' + name + ' mean')
            f.write('; reduced ' + name + ' st.dev.')
        f.write('\n')
        for ii in range(len(bins) - 1):
            f.write('%.3f' % numpy.mean([bins[ii], bins[ii + 1]]))
            for jj in range(len(means)):
                f.write('\t' + '%.3f' % float(means[jj][ii]))
                f.write('\t' + '%.3f' % float(stds[jj][ii]))
            for jj in range(len(means)):
                f.write('\t' + '%.3f' % float(reduced_means[jj][ii]))
                f.write('\t' + '%.3f' % float(reduced_stds[jj][ii]))
            f.write('\n')

    # todo: optionally, finish implementing with matplotlib
    # rc_means = [numpy.mean([bins[ii], bins[ii + 1]]) for ii in range(len(bins) - 1)]
    # for ii in range(len(bins) - 1):
    #     plt.plot(rc_means, means)