"""
Utility functions that don't properly fit anywhere else. These are operations that aren't specific to any particular
interface or script.
"""

import pytraj
import shutil
import fileinput
import re
import os
import sys
import math
import numpy
import pickle
import copy
import subprocess
from statsmodels.tsa import stattools

def check_commit(filename, settings):
    """
    Check commitment of coordinate file to basins defined by settings.commit_fwd and settings.commit_bwd

    Parameters
    ----------
    filename : str
        Name of .rst7-formatted coordinate file to be checked
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    commit_flag : str
        Either 'fwd' or 'bwd' if the coordinates are in the corresponding basin, or '' if in neither

    """

    traj = pytraj.iterload(filename, settings.topology)
    commit_flag = ''    # initialize

    for i in range(len(settings.commit_fwd[2])):
        if settings.commit_fwd[3][i] == 'lt':
            if pytraj.distance(traj, '@' + str(settings.commit_fwd[0][i]) + ' @' + str(settings.commit_fwd[1][i]),
                               n_frames=1)[0] <= settings.commit_fwd[2][i]:
                commit_flag = 'fwd'     # if a committor test is passed, testing moves on to the next one.
            else:
                commit_flag = ''
                break                   # if a committor test is not passed, all testing in this direction fails
        elif settings.commit_fwd[3][i] == 'gt':
            if pytraj.distance(traj, '@' + str(settings.commit_fwd[0][i]) + ' @' + str(settings.commit_fwd[1][i]),
                               n_frames=1)[0] >= settings.commit_fwd[2][i]:
                commit_flag = 'fwd'
            else:
                commit_flag = ''
                break
        else:
            raise ValueError('An incorrect committor definition \"' + settings.commit_fwd[3][i] + '\" was given for '
                             'index ' + str(i) + ' in the \'fwd\' direction.')

    if commit_flag == '':               # only bother checking for bwd commitment if not fwd committed
        for i in range(len(settings.commit_bwd[2])):
            if settings.commit_bwd[3][i] == 'lt':
                if pytraj.distance(traj, '@' + str(settings.commit_bwd[0][i]) + ' @' + str(settings.commit_bwd[1][i]),
                                   n_frames=1)[0] <= settings.commit_bwd[2][i]:
                    commit_flag = 'bwd' # if a committor test is passed, testing moves on to the next one.
                else:
                    commit_flag = ''
                    break               # if a committor test is not passed, all testing in this direction fails
            elif settings.commit_bwd[3][i] == 'gt':
                if pytraj.distance(traj, '@' + str(settings.commit_bwd[0][i]) + ' @' + str(settings.commit_bwd[1][i]),
                                   n_frames=1)[0] >= settings.commit_bwd[2][i]:
                    commit_flag = 'bwd'
                else:
                    commit_flag = ''
                    break
            else:
                raise ValueError('An incorrect committor definition \"' + settings.commit_bwd[3][i] + '\" was given for'
                                 ' index ' + str(i) + ' in the \'bwd\' direction.')

    return commit_flag


def get_cvs(filename, settings, reduce=False):
    """
    Get CV values for a coordinate file given by filename, as well as rates of change if settings.include_qdot = True.

    If reduce = True, the returned CVs will be reduced to between 0 and 1 based on the minimum and maximum values of
    that CV in as.out, which is assumed to exist at settings.as_out_file.
    Parameters
    ----------
    filename : str
        Name of .rst7-formatted coordinate file to be checked
    settings : argparse.Namespace
        Settings namespace object
    reduce : bool
        Boolean determining whether to reduce the CV values for use in evaluating an RC value that uses reduced values

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
                             ' else by that file being truncated. The offending file is: ' + filename)

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
        return (float(unreduced_value) - this_min) / (this_max - this_min)

    traj = pytraj.iterload(filename, settings.topology)

    rc_minmax = [[],[]]
    if reduce:
        # Prepare cv_minmax list
        asout_lines = [[float(item) for item in line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')] for line in open(settings.as_out_file, 'r').readlines()]
        open('as.out', 'r').close()
        mapped = list(map(list, zip(*asout_lines)))
        rc_minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]

    output = ''
    values = []
    local_index = -1
    for cv in settings.cvs:
        local_index += 1
        evaluation = eval(cv)
        if settings.include_qdot:  # want to save values for later
            values.append(float(evaluation))
        if reduce:    # legacy from original atesa, for reducing values to between 0 and 1
            evaluation = reduce_cv(evaluation, local_index, rc_minmax)
        output += str(evaluation) + ' '
    if settings.include_qdot:  # if True, then we want to include rate of change for every CV, too
        # Strategy here is to write a new temporary .rst7 file by incrementing all the coordinate values by their
        # corresponding velocity values, load it as a new iterload object, and then rerun our analysis on that.
        incremented_filename = increment_coords()
        traj = pytraj.iterload(incremented_filename, settings.topology)
        local_index = -1
        for cv in settings.cvs:
            local_index += 1
            evaluation = eval(cv) - values[local_index]  # Subtract value 1/20.455 ps earlier from value of cv
            if reduce:
                evaluation = reduce_cv(evaluation, local_index + len(settings.cvs), rc_minmax)
            output += str(evaluation) + ' '
        os.remove(incremented_filename)     # clean up temporary file

    if output[-1] == ' ':
        output = output[:-1]    # remove trailing space

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
    n_atoms = pattern.findall(byline[1])[0]     # number of atoms indicated on second line of .rst7 file
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
        rc_definition = rc_definition.replace('CV' + str(i), str(cv_list[i]))

    # Evaluate the filled-in rc_definition and return the result
    return eval(rc_definition)


def resample(settings, suffix='', write_raw=True):
    """
    Resample each shooting point in each thread with different CV definitions to produce new as.out files with extant
    aimless shooting data.

    Produces two files: as_raw.out, containing CVs from every shooting point, and as_decorr.out, containing only those
    shooting points where every CV in each thread has decorrelated from its initial coordinates.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object
    suffix : str
        String to append to end of as_decorr.out file name(s) (before extension) to differentiate them
    write_raw : bool
        If True, writes a new copy of as_raw.out; if false, only writes as_decorr.out (if applicable). Also if True and
        settings.information_error_checking = True, writes new partial as_raw.out files and invokes information_error.py
        to build a new info_err.out file.

    Returns
    -------
    None

    """

    # todo: test this thoroughly using a dummy thread and a manual decorrelation time calculation using different software

    # This function is sometimes called from outside the working directory, so make sure we're there
    os.chdir(settings.working_directory)

    # Remove pre-existing as_raw.out if any, initialize new one
    if write_raw:
        open(settings.working_directory + '/as_raw.out', 'w').close()
        if settings.information_error_checking:
            open(settings.working_directory + '/as_raw_timestamped.out', 'w').close()
            open(settings.working_directory + '/info_err.out', 'w').close()
    open(settings.working_directory + '/as_decorr' + suffix + '.out', 'w').close()

    # Load in allthreads from restart.pkl
    try:
        allthreads = pickle.load(open('restart.pkl', 'rb'))
    except FileNotFoundError:
        raise FileNotFoundError('resample = True requires restart.pkl, but could not find one in working directory: '
                                + settings.working_directory)

    # iterate through each thread's history.init_coords list
    for thread in allthreads:
        thread.this_cvs_list = []       # initialize full nested list of CV values for this thread
        thread.cvs_for_later = []       # need this one with empty lists for failed moves, for indexing reasons
        for step_index in range(len(thread.history.prod_results)):
            if thread.history.prod_results[step_index][0] in ['fwd', 'bwd']:
                if thread.history.prod_results[step_index][0] == 'fwd':
                    this_basin = 'B'
                else:  # 'bwd'
                    this_basin = 'A'

                # Get CVs for this shooting point
                try:
                    this_cvs = get_cvs(thread.history.init_coords[step_index][0], settings)

                    # Write CVs to as_raw.out
                    if write_raw:
                        open(settings.working_directory + '/as_raw.out', 'a').write(this_basin + ' <- ')
                        open(settings.working_directory + '/as_raw.out', 'a').write(this_cvs + '\n')
                        open(settings.working_directory + '/as_raw.out', 'a').close()
                        if settings.information_error_checking:
                            open(settings.working_directory + '/as_raw_timestamped.out', 'a').write(str(thread.history.timestamps[step_index]) + ' ')
                            open(settings.working_directory + '/as_raw_timestamped.out', 'a').write(this_basin + ' <- ')
                            open(settings.working_directory + '/as_raw_timestamped.out', 'a').write(this_cvs + '\n')
                            open(settings.working_directory + '/as_raw_timestamped.out', 'a').close()

                    # Append this_cvs to running list for evaluating decorrelation time
                    thread.this_cvs_list.append([[float(item) for item in this_cvs.split(' ')], thread.history.timestamps[step_index]])
                    thread.cvs_for_later.append([float(item) for item in this_cvs.split(' ')])
                except IndexError:  # getting cv's failed (maybe corrupt coordinate file) so consider this step failed
                    thread.cvs_for_later.append([])
            else:
                thread.cvs_for_later.append([])

    if write_raw and settings.information_error_checking:   # sort timestamped output file
        shutil.copy(settings.working_directory + '/as_raw_timestamped.out', settings.working_directory + '/as_raw_timestamped_copy' + suffix + '.out')
        open(settings.working_directory + '/as_raw_timestamped.out', 'w').close()
        with open(settings.working_directory + '/as_raw_timestamped_copy' + suffix + '.out', 'r') as f:
            for line in sorted(f):
                open(settings.working_directory + '/as_raw_timestamped.out', 'a').write(line)
            open(settings.working_directory + '/as_raw_timestamped.out', 'a').close()
        os.remove(settings.working_directory + '/as_raw_timestamped_copy' + suffix + '.out')

    # Construct list of data lengths to perform decorrelation for
    if write_raw and settings.information_error_checking:
        lengths = [leng for leng in range(settings.information_error_freq, len(open(settings.working_directory + '/as_raw_timestamped.out', 'r').readlines()) + 1, settings.information_error_freq)]
        pattern = re.compile('[0-9]+')  # pattern for reading out timestamp from string
    else:
        lengths = [len(open(settings.working_directory + '/as_raw.out', 'r').readlines())]
        pattern = None

    # Assess decorrelation and write as_decorr.out
    for length in lengths:
        if write_raw and settings.information_error_checking:
            suffix = '_' + str(length)     # only use-case with multiple lengths, so this keeps them from stepping on one another's toes
            open(settings.working_directory + '/as_decorr' + suffix + '.out', 'w').close()
            cutoff_timestamp = int(pattern.findall(open(settings.working_directory + '/as_raw_timestamped.out', 'r').readlines()[length])[0])
        else:
            cutoff_timestamp = math.inf
        for thread in allthreads:
            if thread.this_cvs_list:       # if there were any 'fwd' or 'bwd' results in this thread
                mapped = list(map(list, zip(*[item[0] for item in thread.this_cvs_list if item[1] <= cutoff_timestamp])))   # list of lists of values of each CV

                slowest_lag = -1    # initialize running tally of slowest autocorrelation time among CVs in this thread
                if settings.include_qdot:
                    ndims = len(thread.this_cvs_list[0]) / 2   # number of non-rate-of-change CVs
                    if not ndims % 1 == 0:
                        raise ValueError('include_qdot = True but an odd number of dimensions were found in the threads'
                                         ' in restart.pkl')
                    ndims = int(ndims)
                else:
                    ndims = len(thread.this_cvs_list[0])

                for dim_index in range(ndims):
                    this_cv = mapped[dim_index]
                    this_autocorr = [stattools.acf(this_cv[lag:], nlags=1, fft=True)[-1] for lag in range(len(this_cv) - 1)]
                    in_bounds = [this_autocorr[lag] < 1.96 / numpy.sqrt(len(this_cv[lag:])) for lag in range(len(this_autocorr))]
                    running = 0             # running count of successive steps "in bounds"
                    begin = len(in_bounds)  # initialize begin with maximum possible value
                    this_index = -1         # running count of each step regardless of whether "in bounds" or not
                    for item in in_bounds:
                        this_index += 1
                        if (not item and running > 0) or this_index == len(in_bounds) - 1:  # leaves threshold (in -> out) or ends
                            if running >= begin:    # only get slowest_lag if been in bounds for at least as long as it took to get in bounds
                                if begin > slowest_lag:     # only update slowest_lag if it's slower than the previous one
                                    slowest_lag = begin
                                break       # once we've found the lag for this CV move on to the next one
                            begin = 0
                            running = 0
                        elif item and running == 0:  # entered threshold after being outside (out -> in)
                            running += 1
                            begin = this_index  # index where we begin being "in bounds"
                        elif item:  # inside threshold already and still in (in -> in)
                            running += 1
                        else:  # was not in threshold and continues not to be (out -> out)
                            begin = 0
                            running = 0

                if slowest_lag > 0:     # only proceed to writing to as_decorr.out if a valid slowest_lag was found
                    # Write to as_decorr.out the same way as to as_raw.out above, but starting the range at slowest_lag
                    for step_index in range(slowest_lag, len(thread.history.prod_results)):
                        if thread.history.prod_results[step_index][0] in ['fwd', 'bwd']:
                            if thread.history.prod_results[step_index][0] == 'fwd':
                                this_basin = 'B'
                            else:  # 'bwd'
                                this_basin = 'A'

                            # Get CVs for this shooting point
                            if thread.cvs_for_later[step_index]:
                                this_cvs = thread.cvs_for_later[step_index]    # retrieve CVs from last evaluation

                                # Write CVs to as_raw.out
                                open(settings.working_directory + '/as_decorr' + suffix + '.out', 'a').write(this_basin + ' <- ')
                                open(settings.working_directory + '/as_decorr' + suffix + '.out', 'a').write(' '.join([str(item) for item in this_cvs]) + '\n')
                                open(settings.working_directory + '/as_decorr' + suffix + '.out', 'a').close()

        if write_raw and settings.information_error_checking:
            # Add resample_override to settings object to avoid information_error.py calling this function
            settings.resample_override = True
            temp_settings = copy.deepcopy(settings)     # initialize temporary copy of settings to modify
            temp_settings.__dict__.pop('env')           # env attribute is not picklable
            pickle.dump(temp_settings, open('settings.pkl', 'wb'), protocol=2)

            command = 'information_error.py as_decorr' + suffix + '.out'
            process = subprocess.check_call(command.split(' '), stdout=sys.stdout, preexec_fn=os.setsid)
