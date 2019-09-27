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
    that CV in as.out, which is assumed to exist in the present directory. If it does not, get_cvs will raise a
    FileNotFoundError.

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
        asout_lines = [[float(item) for item in line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(' ')] for line in open('as.out', 'r').readlines()]
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


def resample(allthreads, settings):
    """
    Resample each shooting point in each thread with different CV definitions to produce a new as.out file with extant
    aimless shooting data.

    Parameters
    ----------
    allthreads : list
        List of every thread to resample
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    pass    # todo: implement
