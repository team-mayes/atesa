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


def get_cvs(filename, settings):
    """
    Get CV values for a coordinate file given by filename, as well as rates of change if settings.include_qdot = True

    Parameters
    ----------
    filename : str
        Name of .rst7-formatted coordinate file to be checked
    settings : argparse.Namespace
        Settings namespace object

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

    traj = pytraj.iterload(filename, settings.topology)

    output = ''
    values = []
    local_index = -1
    for cv in settings.cvs:
        local_index += 1
        evaluation = eval(cv)
        if settings.include_qdot:  # want to save values for later
            values.append(float(evaluation))
        # if reduce:    # legacy from original atesa, for reducing values to between 0 and 1
        #     evaluation = reduce_cv(evaluation, local_index)
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
            # if reduce:
            #     evaluation = reduce_cv(evaluation, local_index + len(settings.candidateops))
            output += str(evaluation) + ' '
        os.remove(incremented_filename)     # clean up temporary file

    return output
